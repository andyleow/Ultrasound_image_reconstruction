%%Tested platform

%Windows 
Matlab 2018b
CUDA 8.0
Visual Studio 2015

matlab compile commmand

If (CUDA toolkits 9.1 installed) 
    mexcuda -v cuDAS_mainv1.cu cuDAS_functionV1.cu
else
    (windows)
    mex -largeArrayDims -f nvcc_msvcpp2015.xml NVCC_FLAGS="" -v cuDAS_mainv1.cu cuDAS_functionV1.cu	

%Ubuntu system
Matlab 2017b
CUDA 8.0
   mex -largeArrayDims -f nvcc9_g++.xml NVCC_FLAGS="" -v cuDAS_int.cu cuDAS_function.cu 
CUDA 9.0
   mex -largeArrayDims -f nvcc9_g++.xml NVCC_FLAGS="" -v cuDAS_int.cu cuDAS_function.cu 	
	

** Parrallel computing toolbox must be installed on matlab in order to use mexcuda function

# Possible error 
g++: error: /usr/local/MATLAB/R2017b/bin/glnxa64/libcudart.so.9.0: No such file or directory
Fix: Copy "libcudart.so.9.0" and "libcufft.so.9.0" from "/usr/local/cuda/lib64" and paste to "/usr/local/MATLAB/R2017b/bin/glnxa64" by entering the command below:
sudo -i 
cp /usr/local/cuda/lib64/libcudart.so.9.0 /usr/local/MATLAB/R2017b/bin/glnxa64
cp /usr/local/cuda/lib64/libcufft.so.9.0 /usr/local/MATLAB/R2017b/bin/glnxa64

#In ubuntu, Matlab compiler makeFile, the following need to be set: 
DEFINES="--compiler-options=-D_GNU_SOURCE,$MATLABMEX -DCUDA_API_PER_THREAD_DEFAULT_STREAM  -use_fast_math"
NVCCFLAGS="-gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=&#92;&quot;sm_61,compute_61&#92;&quot; -std=c++11 $NVCC_FLAGS" 
CXXOPTIMFLAGS="-O2 -DNDEBUG" 
LDOPTIMFLAGS="-O2"

