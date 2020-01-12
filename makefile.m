%%Function to Compiling source code
addpath(genpath('./src'));

if isunix
    try
        %Matlab legacy compiler
        mexcuda -Llibcufft.so -v src/gpuDAS/cuDAS_int.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_int
        mexcuda -Llibcufft.so -v src/gpuDAS/cuDAS_single.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_single
        mexcuda -Llibcufft.so -v src/gpuFBF/cuFBF_int.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_int
        mexcuda -Llibcufft.so -v src/gpuFBF/cuFBF_single.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_single
    catch
        %Custom compiler
        mex -largeArrayDims -f src/gpuDAS/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_int.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_int
        mex -largeArrayDims -f src/gpuDAS/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_single.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_single
        mex -largeArrayDims -f src/gpuFBF/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_int.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_int
        mex -largeArrayDims -f src/gpuFBF/nvcc9_g++.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_single.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_singl
    end
elseif ispc
    try
        %Matlab legacy compiler
        mexcuda -Lcufft.lib src/gpuDAS/nvcc_msvcpp2015.xml -v src/gpuDAS/cuDAS_int.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_int
        mexcuda -Lcufft.lib src/gpuDAS/nvcc_msvcpp2015.xml -v src/gpuDAS/cuDAS_single.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_single
        mexcuda -Lcufft.lib src/gpuFBF/nvcc_msvcpp2015.xml -v src/gpuFBF/cuFBF_int.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_int
        mexcuda -Lcufft.lib src/gpuFBF/nvcc_msvcpp2015.xml -v src/gpuFBF/cuFBF_single.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_single
        
    catch
        %Custom compiler
        mex -largeArrayDims -f src/gpuDAS/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_int.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_int
        mex -largeArrayDims -f src/gpuDAS/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuDAS/cuDAS_single.cu src/gpuDAS/cuDAS_function.cu -output src/gpuDAS/cuDAS_single
        mex -largeArrayDims -f src/gpuFBF/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_int.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_int
        mex -largeArrayDims -f src/gpuFBF/nvcc_msvcpp2015.xml NVCC_FLAGS="" -v src/gpuFBF/cuFBF_single.cu src/gpuFBF/cuFBF_function.cu -output src/gpuFBF/cuFBF_single
    end
else
    error("Compiling source code failed. Non linux or windows operating system");
end
