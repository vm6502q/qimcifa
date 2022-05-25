//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2022. All rights reserved.
//
// This is a multithreaded, universal quantum register simulation, allowing
// (nonphysical) register cloning and direct measurement of probability and
// phase, to leverage what advantages classical emulation of qubits can have.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#pragma once

#define _USE_MATH_DEFINES

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#endif

#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#if defined(__APPLE__)
#define CL_SILENCE_DEPRECATION
#include <OpenCL/cl.hpp>
#elif defined(_WIN32) || ENABLE_SNUCL
#include <CL/cl.hpp>
#else
#include <CL/cl2.hpp>
#endif

namespace Qimcifa {

class OCLDeviceCall;

class OCLDeviceContext;

typedef std::shared_ptr<OCLDeviceContext> DeviceContextPtr;
typedef std::shared_ptr<std::vector<cl::Event>> EventVecPtr;

enum OCLAPI {
    OCL_API_UNKNOWN = 0,
    OCL_API_QIMCIFA_BATCH
};

struct OCLKernelHandle {
    OCLAPI oclapi;
    std::string kernelname;

    OCLKernelHandle(OCLAPI o, std::string kn)
        : oclapi(o)
        , kernelname(kn)
    {
    }
};

class OCLDeviceCall {
protected:
    std::lock_guard<std::mutex> guard;

public:
    // A cl::Kernel is unique object which should always be taken by reference, or the OCLDeviceContext will lose
    // ownership.
    cl::Kernel& call;
    OCLDeviceCall(const OCLDeviceCall&);

protected:
    OCLDeviceCall(std::mutex& m, cl::Kernel& c)
        : guard(m)
        , call(c)
    {
    }

    friend class OCLDeviceContext;

private:
    OCLDeviceCall& operator=(const OCLDeviceCall&) = delete;
};

class OCLDeviceContext {
public:
    cl::Platform platform;
    cl::Device device;
    cl::Context context;
    int64_t context_id;
    int64_t device_id;
    cl::CommandQueue queue;
    EventVecPtr wait_events;

protected:
    std::mutex waitEventsMutex;
    std::map<OCLAPI, cl::Kernel> calls;
    std::map<OCLAPI, std::unique_ptr<std::mutex>> mutexes;

private:
    const size_t procElemCount = device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
    const size_t maxWorkItems = device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0];
    const size_t maxAlloc = device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();
    const size_t globalSize = device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
    size_t preferredSizeMultiple;
    size_t preferredConcurrency;

public:
    OCLDeviceContext(cl::Platform& p, cl::Device& d, cl::Context& c, int64_t dev_id, int64_t cntxt_id)
        : platform(p)
        , device(d)
        , context(c)
        , context_id(cntxt_id)
        , device_id(dev_id)
        , preferredSizeMultiple(0U)
        , preferredConcurrency(0U)
    {
        cl_int error;
        queue = cl::CommandQueue(context, d, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &error);
        if (error != CL_SUCCESS) {
            queue = cl::CommandQueue(context, d);
        }

        wait_events =
            std::shared_ptr<std::vector<cl::Event>>(new std::vector<cl::Event>(), [](std::vector<cl::Event>* vec) {
                vec->clear();
                delete vec;
            });
    }

    OCLDeviceCall Reserve(OCLAPI call) { return OCLDeviceCall(*(mutexes[call]), calls[call]); }

    EventVecPtr ResetWaitEvents()
    {
        std::lock_guard<std::mutex> guard(waitEventsMutex);
        EventVecPtr waitVec = std::move(wait_events);
        wait_events =
            std::shared_ptr<std::vector<cl::Event>>(new std::vector<cl::Event>(), [](std::vector<cl::Event>* vec) {
                vec->clear();
                delete vec;
            });
        return waitVec;
    }

    void LockWaitEvents() { waitEventsMutex.lock(); }

    void UnlockWaitEvents() { waitEventsMutex.unlock(); }

    void WaitOnAllEvents()
    {
        std::lock_guard<std::mutex> guard(waitEventsMutex);
        if ((wait_events.get())->size()) {
            cl::Event::waitForEvents((const std::vector<cl::Event>&)*(wait_events.get()));
            wait_events->clear();
        }
    }

    size_t GetPreferredSizeMultiple()
    {
        return preferredSizeMultiple
            ? preferredSizeMultiple
            : preferredSizeMultiple =
                  calls[OCL_API_QIMCIFA_BATCH].getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(
                      device);
    }

    size_t GetPreferredConcurrency()
    {
        if (preferredConcurrency) {
            return preferredConcurrency;
        }

        preferredConcurrency = procElemCount * GetPreferredSizeMultiple();

        return preferredConcurrency;
    }

    size_t GetProcElementCount() { return procElemCount; }
    size_t GetMaxWorkItems() { return maxWorkItems; }
    size_t GetMaxAlloc() { return maxAlloc; }
    size_t GetGlobalSize() { return globalSize; }

    friend class OCLEngine;
};

struct InitOClResult {
    std::vector<DeviceContextPtr> all_dev_contexts;
    DeviceContextPtr default_dev_context;

    InitOClResult()
        : all_dev_contexts()
        , default_dev_context(NULL)
    {
        // Intentionally left blank
    }

    InitOClResult(std::vector<DeviceContextPtr> adc, DeviceContextPtr ddc)
        : all_dev_contexts(adc)
        , default_dev_context(ddc)
    {
        // Intentionally left blank
    }
};

/** "Qimcifa::OCLEngine" manages the single OpenCL context. */
class OCLEngine {
public:
    // See https://stackoverflow.com/questions/1008019/c-singleton-design-pattern
    /// Get a pointer to the Instance of the singleton. (The instance will be instantiated, if it does not exist yet.)
    static OCLEngine& Instance()
    {
        static OCLEngine instance;
        return instance;
    }
    /// Get default location for precompiled binaries:
    static std::string GetDefaultBinaryPath()
    {
#if ENABLE_ENV_VARS
        if (getenv("QIMCIFA_OCL_PATH")) {
            std::string toRet = std::string(getenv("QIMCIFA_OCL_PATH"));
            if ((toRet.back() != '/') && (toRet.back() != '\\')) {
#if defined(_WIN32) && !defined(__CYGWIN__)
                toRet += "\\";
#else
                toRet += "/";
#endif
            }
            return toRet;
        }
#endif
#if defined(_WIN32) && !defined(__CYGWIN__)
        return std::string(getenv("HOMEDRIVE") ? getenv("HOMEDRIVE") : "") +
            std::string(getenv("HOMEPATH") ? getenv("HOMEPATH") : "") + "\\.qimcifa\\";
#else
        return std::string(getenv("HOME") ? getenv("HOME") : "") + "/.qimcifa/";
#endif
    }
    /// Initialize the OCL environment, with the option to save the generated binaries. Binaries will be saved/loaded
    /// from the folder path "home". This returns a Qimcifa::OCLInitResult object which should be passed to
    /// SetDeviceContextPtrVector().
    static InitOClResult InitOCL(bool buildFromSource = false, bool saveBinaries = false, std::string home = "*");

    /// Get a pointer one of the available OpenCL contexts, by its index in the list of all contexts.
    DeviceContextPtr GetDeviceContextPtr(const int64_t& dev = -1);
    /// Get the list of all available devices (and their supporting objects).
    std::vector<DeviceContextPtr> GetDeviceContextPtrVector();
    /** Set the list of DeviceContextPtr object available for use. If one takes the result of
     * GetDeviceContextPtrVector(), trims items from it, and sets it with this method, (at initialization, before any
     * QEngine objects depend on them,) all resources associated with the removed items are freed.
     */
    void SetDeviceContextPtrVector(std::vector<DeviceContextPtr> vec, DeviceContextPtr dcp = nullptr);
    /// Get the count of devices in the current list.
    int GetDeviceCount() { return all_device_contexts.size(); }
    /// Get default device ID.
    size_t GetDefaultDeviceID() { return default_device_context->device_id; }
    /// Pick a default device, for QEngineOCL instances that don't specify a preferred device.
    void SetDefaultDeviceContext(DeviceContextPtr dcp);

    OCLEngine(OCLEngine const&) = delete;
    void operator=(OCLEngine const&) = delete;

private:
    static const std::vector<OCLKernelHandle> kernelHandles;
    static const std::string binary_file_prefix;
    static const std::string binary_file_ext;

    std::vector<DeviceContextPtr> all_device_contexts;
    DeviceContextPtr default_device_context;

    OCLEngine(); // Private so that it can  not be called

    /// Make the program, from either source or binary
    static cl::Program MakeProgram(bool buildFromSource, std::string path, std::shared_ptr<OCLDeviceContext> devCntxt);
    /// Save the program binary:
    static void SaveBinary(cl::Program program, std::string path, std::string fileName);
};

} // namespace Qimcifa
