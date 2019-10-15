#  Real time GPU Beamformer for ultrafast ultrasound imaging 

This repository contains GPU-accelerated codes for real-time reconstruction of plane wave images. 
These includes: 
	*CPU delay-and-sum beamformer
	*CPU Fourier beamformer
	*GPU delay-and-sum beamformer
	*GPU Fourier beamformer

The GPU programs are written in C/CUDA but interface with matlab using mex function for real-time visualisation and integration with Vantage research ultrasound scanner

# Requirements
The program works for both Windows and linux system and has been tested using:
* matlab 2017b 
* cuda toolkits 9.0
* visual studio 2015 (windows)

#  Installation 
Clone this repository 

```
git clone https://github.com/andyleow/Ultrasound_image_reconstruction.git
```

To compile the GPU from source code, user can either execute the "makefile.m" from matlab or independently execute the following commands in the src file directory 
in Matlab
* Windows 
```
mex -largeArrayDims -f nvcc_msvcpp2015.xml NVCC_FLAGS="" -v cuDAS_int.cu cuDAS_function.cu
mex -largeArrayDims -f nvcc_msvcpp2015.xml NVCC_FLAGS="" -v cuDAS_single.cu cuDAS_function.cu
mex -largeArrayDims -f nvcc_msvcpp2015.xml NVCC_FLAGS="" -v cuFBF_int.cu cuFBF_functionV1.cu
mex -largeArrayDims -f nvcc_msvcpp2015.xml NVCC_FLAGS="" -v cuFBF_single.cu cuFBF_functionV1.cu

```

* Linux 
```
mex -largeArrayDims -f nvcc9_g++.xml NVCC_FLAGS="" -v cuDAS_int.cu cuDAS_function.cu 
mex -largeArrayDims -f nvcc9_g++.xml NVCC_FLAGS="" -v cuDAS_single.cu cuDAS_function.cu 
mex -largeArrayDims -f nvcc9_g++.xml NVCC_FLAGS="" -v cuFBF_int.cu cuFBF_function.cu 
mex -largeArrayDims -f nvcc9_g++.xml NVCC_FLAGS="" -v cuFBF_single.cu cuFBF_function.cu 
```

# Demonstration
Offline example can be found in`Example/Offline` 
Online execution on Vantage scanner can be found on `Example/Vantage`

# Citation 


# Future work




