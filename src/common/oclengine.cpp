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

#include "common/oclengine.hpp"

#include "kiss09cl.hpp"
#include "big_integercl.hpp"
#include "qimcifa_uint64cl.hpp"

#include <algorithm>
#include <iostream>
#include <memory>

namespace Qimcifa {

/// "Qimcifa::OCLEngine" manages the single OpenCL context

// Public singleton methods to get pointers to various methods
DeviceContextPtr OCLEngine::GetDeviceContextPtr(const int64_t& dev)
{
    if ((dev >= GetDeviceCount()) || (dev < -1)) {
        throw "Invalid OpenCL device selection";
    } else if (dev == -1) {
        return default_device_context;
    } else {
        return all_device_contexts[dev];
    }
}

// clang-format off
const std::vector<OCLKernelHandle> OCLEngine::kernelHandles = {
    OCLKernelHandle(OCL_API_QIMCIFA_BATCH, "qimcifa_batch")
};
// clang-format on

const std::string OCLEngine::binary_file_prefix("qimcifa_ocl_dev_");
const std::string OCLEngine::binary_file_ext(".ir");

std::vector<DeviceContextPtr> OCLEngine::GetDeviceContextPtrVector() { return all_device_contexts; }
void OCLEngine::SetDeviceContextPtrVector(std::vector<DeviceContextPtr> vec, DeviceContextPtr dcp)
{
    all_device_contexts = vec;
    if (dcp != nullptr) {
        default_device_context = dcp;
    }
}

void OCLEngine::SetDefaultDeviceContext(DeviceContextPtr dcp) { default_device_context = dcp; }

cl::Program OCLEngine::MakeProgram(bool buildFromSource, std::string path, std::shared_ptr<OCLDeviceContext> devCntxt)
{
    FILE* clBinFile;
    cl::Program program;
    cl_int buildError = -1;
    std::vector<int> binaryStatus;
    if (!buildFromSource && (clBinFile = fopen(path.c_str(), "r"))) {
        struct stat statSize;
        if (fstat(fileno(clBinFile), &statSize)) {
            std::cout << "Binary error: Invalid file fstat result. (Falling back to JIT.)" << std::endl;
        } else {
            size_t lSize = statSize.st_size;
            std::vector<unsigned char> buffer(lSize);
            size_t lSizeResult = fread(&buffer[0U], sizeof(unsigned char), lSize, clBinFile);
            fclose(clBinFile);

            if (lSizeResult != lSize) {
                std::cout << "Binary warning: Binary file size and read result length do not match. (Attempting to "
                             "build anyway.)"
                          << std::endl;
            }

#if defined(__APPLE__) || (defined(_WIN32) && !defined(__CYGWIN__)) || ENABLE_SNUCL
            program = cl::Program(devCntxt->context, { devCntxt->device },
                { std::pair<const void*, size_t>(&buffer[0U], buffer.size()) }, &binaryStatus, &buildError);
#else
            program = cl::Program(devCntxt->context, { devCntxt->device }, { buffer }, &binaryStatus, &buildError);
#endif

            if ((buildError != CL_SUCCESS) || (binaryStatus[0U] != CL_SUCCESS)) {
                std::cout << "Binary error: " << buildError << ", " << binaryStatus[0U] << " (Falling back to JIT.)"
                          << std::endl;
            } else {
                std::cout << "Loaded binary from: " << path << std::endl;
            }
        }
    }

    // If, either, there are no cached binaries, or binary loading failed, then fall back to JIT.
    if (buildError == CL_SUCCESS) {
        return program;
    }
    
    // TODO: This needs manual file I/O.
    // create the programs that we want to execute on the devices
    cl::Program::Sources sources;
    sources.push_back({ (const char*)kiss09_cl, (long unsigned int)kiss09_cl_len });
    sources.push_back({ (const char*)big_integer_cl, (long unsigned int)big_integer_cl_len });
    sources.push_back({ (const char*)qimcifa_uint64_cl, (long unsigned int)qimcifa_uint64_cl_len });

    program = cl::Program(devCntxt->context, sources);
    std::cout << "Built JIT." << std::endl;
    
    return program;
}

void OCLEngine::SaveBinary(cl::Program program, std::string path, std::string fileName)
{
    std::vector<size_t> clBinSizes = program.getInfo<CL_PROGRAM_BINARY_SIZES>();
    size_t clBinSize = 0U;
    int64_t clBinIndex = 0;

    for (size_t i = 0U; i < clBinSizes.size(); i++) {
        if (clBinSizes[i]) {
            clBinSize = clBinSizes[i];
            clBinIndex = i;
            break;
        }
    }

    std::cout << "Binary size:" << clBinSize << std::endl;

#if defined(_WIN32) && !defined(__CYGWIN__)
    int err = _mkdir(path.c_str());
#else
    int err = mkdir(path.c_str(), 0700);
#endif
    if (err != -1) {
        std::cout << "Making directory: " << path << std::endl;
    }

    FILE* clBinFile = fopen((path + fileName).c_str(), "w");
#if defined(__APPLE__) || (defined(_WIN32) && !defined(__CYGWIN__)) || ENABLE_SNUCL
    std::vector<char*> clBinaries = program.getInfo<CL_PROGRAM_BINARIES>();
    char* clBinary = clBinaries[clBinIndex];
    fwrite(clBinary, clBinSize, sizeof(char), clBinFile);
#else
    std::vector<std::vector<unsigned char>> clBinaries = program.getInfo<CL_PROGRAM_BINARIES>();
    std::vector<unsigned char> clBinary = clBinaries[clBinIndex];
    fwrite(&clBinary[0U], clBinSize, sizeof(unsigned char), clBinFile);
#endif
    fclose(clBinFile);
}

InitOClResult OCLEngine::InitOCL(bool buildFromSource, bool saveBinaries, std::string home)
{

    if (home == "*") {
        home = GetDefaultBinaryPath();
    }
    // get all platforms (drivers), e.g. NVIDIA

    std::vector<cl::Platform> all_platforms;
    std::vector<cl::Device> all_devices;
    std::vector<int64_t> device_platform_id;
    cl::Platform default_platform;
    cl::Device default_device;
    std::vector<DeviceContextPtr> all_dev_contexts;
    DeviceContextPtr default_dev_context;

    cl::Platform::get(&all_platforms);

    if (!all_platforms.size()) {
        std::cout << " No platforms found. Check OpenCL installation!\n";
        return InitOClResult();
    }

    // get all devices
    std::vector<cl::Platform> devPlatVec;
    std::vector<std::vector<cl::Device>> all_platforms_devices;
    for (size_t i = 0U; i < all_platforms.size(); i++) {
        all_platforms_devices.push_back(std::vector<cl::Device>());
        all_platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &(all_platforms_devices[i]));
        for (size_t j = 0U; j < all_platforms_devices[i].size(); j++) {
            // VirtualCL seems to break if the assignment constructor of cl::Platform is used here from the original
            // list. Assigning the object from a new query is always fine, though. (They carry the same underlying
            // platform IDs.)
            std::vector<cl::Platform> temp_platforms;
            cl::Platform::get(&temp_platforms);
            devPlatVec.push_back(temp_platforms[i]);
            device_platform_id.push_back(i);
        }
        all_devices.insert(all_devices.end(), all_platforms_devices[i].begin(), all_platforms_devices[i].end());
    }
    if (!all_devices.size()) {
        std::cout << " No devices found. Check OpenCL installation!\n";
        return InitOClResult();
    }

    int64_t deviceCount = all_devices.size();
    // prefer the last device because that's usually a GPU or accelerator; device[0U] is usually the CPU
    int64_t dev = deviceCount - 1;
    if (getenv("QIMCIFA_OCL_DEFAULT_DEVICE")) {
        dev = std::stoi(std::string(getenv("QIMCIFA_OCL_DEFAULT_DEVICE")));
        if ((dev < 0) || (dev > (deviceCount - 1))) {
            std::cout << "WARNING: Invalid QIMCIFA_OCL_DEFAULT_DEVICE selection. (Falling back to highest index device "
                         "as default.)"
                      << std::endl;
            dev = deviceCount - 1;
        }
    }

    int64_t plat_id = -1;
    std::vector<cl::Context> all_contexts;
    std::vector<std::string> all_filenames;
    for (int64_t i = 0; i < deviceCount; i++) {
        // a context is like a "runtime link" to the device and platform;
        // i.e. communication is possible
        if (device_platform_id[i] != plat_id) {
            plat_id = device_platform_id[i];
            all_contexts.push_back(cl::Context(all_platforms_devices[plat_id]));
        }
        std::shared_ptr<OCLDeviceContext> devCntxt = std::make_shared<OCLDeviceContext>(
            devPlatVec[i], all_devices[i], all_contexts[all_contexts.size() - 1U], i, plat_id);

        std::string fileName = binary_file_prefix + all_devices[i].getInfo<CL_DEVICE_NAME>() + binary_file_ext;
        std::replace(fileName.begin(), fileName.end(), ' ', '_');
        std::string clBinName = home + fileName;

        std::cout << "Device #" << i << ", ";
        cl::Program program = MakeProgram(buildFromSource, clBinName, devCntxt);

        cl_int buildError =
            program.build({ all_devices[i] }, "-cl-strict-aliasing -cl-denorms-are-zero -cl-fast-relaxed-math");
        if (buildError != CL_SUCCESS) {
            std::cout << "Error building for device #" << i << ": " << buildError << ", "
                      << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(all_devices[i])
                      << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(all_devices[i]) << std::endl;

            // The default device was set above to be the last device in the list. If we can't compile for it, we
            // use the first device. If the default is the first device, and we can't compile for it, then we don't
            // have any devices that can compile at all, and the environment needs to be fixed by the user.
            if (i == dev) {
                default_dev_context = all_dev_contexts[0U];
                default_platform = all_platforms[0U];
                default_device = all_devices[0U];
            }

            continue;
        }

        all_dev_contexts.push_back(devCntxt);

        for (unsigned int j = 0U; j < kernelHandles.size(); j++) {
            all_dev_contexts[i]->calls[kernelHandles[j].oclapi] =
                cl::Kernel(program, kernelHandles[j].kernelname.c_str());
            all_dev_contexts[i]->mutexes.emplace(kernelHandles[j].oclapi, new std::mutex);
        }

        std::vector<std::string>::iterator fileNameIt = std::find(all_filenames.begin(), all_filenames.end(), fileName);
        if (saveBinaries && (fileNameIt == all_filenames.end())) {
            std::cout << "OpenCL program #" << i << ", ";
            SaveBinary(program, home, fileName);
        }

        if (i == dev) {
            default_dev_context = all_dev_contexts[i];
            default_platform = all_platforms[plat_id];
            default_device = all_devices[i];
        }
    }

    // For VirtualCL support, the device info can only be accessed AFTER all contexts are created.
    std::cout << "Default platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";
    std::cout << "Default device: #" << dev << ", " << default_device.getInfo<CL_DEVICE_NAME>() << "\n";
    for (int64_t i = 0; i < deviceCount; i++) {
        std::cout << "OpenCL device #" << i << ": " << all_devices[i].getInfo<CL_DEVICE_NAME>() << "\n";
    }

    return InitOClResult(all_dev_contexts, default_dev_context);
}

OCLEngine::OCLEngine()
{
    InitOClResult initResult = InitOCL(false);
    SetDeviceContextPtrVector(initResult.all_dev_contexts, initResult.default_dev_context);
}

} // namespace Qimcifa
