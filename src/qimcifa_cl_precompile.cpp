//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qimcifa contributors 2017-2021. All rights reserved.
//
// This utility precompiles all the OpenCL kernels that Qimcifa needs.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "common/oclengine.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    // Precompile, if OpenCL is available.
    std::cout << "Precompiling OCL kernels..." << std::endl;
    if (argc < 2) {
        std::cout << "Will save to: " << Qimcifa::OCLEngine::GetDefaultBinaryPath() << std::endl;
        Qimcifa::OCLEngine::InitOCL(true, true, Qimcifa::OCLEngine::GetDefaultBinaryPath());
    } else {
        std::cout << "Will save to: " << std::string(argv[1]) << std::endl;
        Qimcifa::OCLEngine::InitOCL(true, true, std::string(argv[1]) + "/");
    }
    std::cout << "Done precompiling OCL kernels." << std::endl;
}
