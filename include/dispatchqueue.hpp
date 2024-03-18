//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// This is a multithreaded, universal quantum register simulation, allowing
// (nonphysical) register cloning and direct measurement of probability and
// phase, to leverage what advantages classical emulation of qubits can have.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

// From https://github.com/embeddedartistry/embedded-resources/blob/master/examples/cpp/dispatch.cpp

#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>

typedef std::function<bool(void)> DispatchFn;

class DispatchQueue {
public:
    DispatchQueue(size_t n)
        : threads_(n)
        , quit_(false)
        , isFinished_(true)
        , isStarted_(false)
        , result(false)
    {
        // Intentionally left blank.
    }
    ~DispatchQueue();

    // Reset output result.
    void resetResult() { result = false; }
    // dispatch and copy
    void dispatch(const DispatchFn& op);
    // finish queue
    bool finish();
    // dump queue
    void dump();
    // check if queue is finished
    bool isFinished() { return isFinished_; }

    // Deleted operations
    DispatchQueue(const DispatchQueue& rhs) = delete;
    DispatchQueue& operator=(const DispatchQueue& rhs) = delete;
    DispatchQueue(DispatchQueue&& rhs) = delete;
    DispatchQueue& operator=(DispatchQueue&& rhs) = delete;

private:
    std::mutex lock_;
    std::vector<std::future<void>> threads_;
    std::queue<DispatchFn> q_;
    std::condition_variable cv_;
    std::condition_variable cvFinished_;
    bool quit_;
    bool isFinished_;
    bool isStarted_;
    bool result;

    void dispatch_thread_handler(void);
};
