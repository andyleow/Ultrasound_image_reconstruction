#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include "mat.h"
#include "mex.h"

// Includes CUDA
#include <cuda_runtime.h>
#include <cuda.h>

typedef signed short int16_t;

//Define total number of stream
#define nStream 4

#define CUDA_ERROR_CHECK
#define cudaCheckErr() __cudaCheckErr( __FILE__, __LINE__)
inline void __cudaCheckErr(char *file, int line, bool abort = true)
{
#ifdef CUDA_ERROR_CHECK
	cudaDeviceSynchronize();
	cudaError err = cudaGetLastError();
	if (err != cudaSuccess)
	{
		mexPrintf("Error");
		fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
		//system("pause");
		if (abort) exit(err);
	}
#endif
}

// Complex data type
typedef float2 Complex;
static __device__ __host__ inline Complex ComplexAdd(Complex, Complex);
static __device__ __host__ inline Complex ComplexScale(Complex, float);
static __device__ __host__ inline Complex ComplexMul(Complex, Complex);
static __global__ void ComplexPointwiseMulAndScale(Complex *, const Complex *, int, int, float);

//Host memory pointer
struct hostPtr {
	int init = 0;
	int numSample;
	int filterLength;
	int padSample;
	int numChannel;
	int numAngle;
	int numFrame;
	int imageHeight;
	int imageWidth;
	int reconMode;
	int gpuID=0;

	float c;
	float fs;
	float ftx;
	float focus;
	float pitch;
	float senseCutoff;

	float *ElePos;
	float *angle;
	float *pixelMapX;
	float *pixelMapZ;
	float *t0;
	float *filter;
	float *LRI_r;
	float *LRI_i;
};

//Cuda memory pointer
struct devicePtr {
	//Pinned memory 
	Complex *h_RFIn;
	Complex *h_LRI;

	//Device memory
	Complex *d_RFIn;
	Complex *d_RFIn_temp;
	Complex *d_LRI;
	Complex *d_filter;
	float *d_pixelMapX;
	float *d_pixelMapZ;
	int *d_coef;

	cudaStream_t *stream;
	cudaArray* d_RFArrayIn;
	cudaMemcpy3DParms myparms = { 0 };

	cufftHandle planFFT_time;
	cufftHandle planFFT_filter;
};

void init_hostPtr(struct hostPtr &hpp, double *h_angle, double *h_ElePos, double *h_filter, double *h_pixelMapX, double *h_pixelMapZ, double *h_t0);

void init_cudaMemory(struct devicePtr &dpp, struct hostPtr hpp);

void DAS_beamforming_GPU_LRI(struct devicePtr &dpp, struct hostPtr hpp);

void DAS_beamforming_GPU_HRI(struct devicePtr &dpp, struct hostPtr hpp);

void DAS_beamforming_GPU_CRI(struct devicePtr &dpp, struct hostPtr hpp);

void free_hostMemory(struct hostPtr &hpp);
void free_cudaMemory(struct devicePtr &dpp);

void firHilbert(Complex *d_RFIn, Complex *d_RFIn_temp, Complex *d_filter, int *d_coef, int numSample, int numChannel, cufftHandle planFFT_time, cudaStream_t stream);

__global__ void EleWiseProduct(float2 *d_RFIn_r, int *d_coef, unsigned int numSample, unsigned int batch);

__global__ void kernel_DAS_PW(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX, const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel, const float t0, const unsigned int sID);

__global__ void kernel_DAS_Diverge(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX, const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel,
	const float t0, const float focus, const unsigned int sID);

__global__ void kernel_DAS_PW_Channel(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX, const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel, const float t0, const unsigned int sID);

__global__ void kernel_DAS_Diverge_Channel(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX , const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel,
	const float t0, const float focus, const unsigned int sID);
