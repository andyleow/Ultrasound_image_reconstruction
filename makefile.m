%%Function to Compiling source code 
addpath(genpath('./src'));

if isunix
    mex -largeArrayDims -f src/gpuDAS/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_int.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_int
    mex -largeArrayDims -f src/gpuDAS/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_single.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_single
    mex -largeArrayDims -f src/gpuFBF/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_int.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_int
    mex -largeArrayDims -f src/gpuFBF/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_single.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_single    
elseif ispc
    mex -largeArrayDims -f src/gpuDAS/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_int.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_int
    mex -largeArrayDims -f src/gpuDAS/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_int.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_single
    mex -largeArrayDims -f src/gpuFBF/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_single.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_int
    mex -largeArrayDims -f src/gpuFBF/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_single.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_single    
else
    error("Compiling source code failed. Non linux or windows operating system");
end
