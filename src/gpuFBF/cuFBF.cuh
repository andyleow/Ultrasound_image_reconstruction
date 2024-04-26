#include <stdio.h>
//#include < stdint.h>
#include <math.h>
#include <cufft.h>
#include "mat.h"
#include "mex.h"

// Includes CUDA
#include <cuda_runtime.h>
#include <cuda.h>

typedef signed short int16_t;

//Define total number of stream
#define nStream 8

#define CUDA_ERROR_CHECK
#define cudaCheckErr() __cudaCheckErr( __FILE__, __LINE__)
inline void __cudaCheckErr(char *file, int line, bool abort = true)
{
#ifdef CUDA_ERROR_CHECK
	cudaDeviceSynchronize();
	cudaError err = cudaGetLastError();
	if (err != cudaSuccess)
	{
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

__global__ void EleWiseProduct(Complex *d_RFIn_r, int *d_coef, unsigned int numSample, unsigned int batch);
__global__ void kernel_complexAxialShift(Complex *d_RFIn, float angle, float t0, unsigned int numSample, int numChannel);
__global__ void kernel_complexLateralShift(Complex *d_RFIn, float *d_kx, unsigned int numSample, int numChannel);

__global__ void kernel_transpose(Complex *odata, Complex *idata, int width, int height);

__global__ void kernel_FourierInterpolation(Complex* __restrict__  d_LRI,
	float* __restrict__ kx, const float angle, const unsigned int numSample, const unsigned int numChannel, const unsigned int sID);


struct hostPtr {
	int init = 0;
	int numSample;
	int filterLength;
	int padLength;
	int padWidth;
	int numChannel;
	int numAngle;
	int numFrame;
	int imageHeight;
	int imageWidth;
	int reconMode;
	int gpuID = 0;

	float c;
	float fs;
	float ftx;
	float focus;
	float pitch;
	float deltaX = 0.f;
	int deltaZ = 0;

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
	Complex *h_RFIn_f;
	Complex *h_LRI_f;

	//Device memory
	Complex *d_RFIn;
	Complex *d_RFIn_tranpose;
	Complex *d_LRI;
	Complex *d_LRI_transpose;
	Complex *d_filter;

	//float *d_angle;
	//float *d_t0;
	int *d_coef;
	float* d_kx;
	cudaStream_t *stream;
	cudaArray* d_RFArrayIn;
	cudaMemcpy3DParms myparms = { 0 };

	cufftHandle planFFT_time_fwd;
	cufftHandle planFFT_time_rvs;
	cufftHandle planFFT_space_fwd;
	cufftHandle planFFT_space_rvs;
	cufftHandle planFFT_filter;

};

void init_hostPtr(struct hostPtr &hpp, double *h_angle, double *h_ElePos, double *h_filter, double *h_pixelMapX, double *h_pixelMapZ, double *h_t0);

void init_cudaMemory(struct devicePtr &dpp, struct hostPtr hpp);

//struct devicePtr init_cudaMemory(int mode, const float *h_filter, float h_pitch, float h_c, float h_fs, int const numFrame,
//	int const numAngle, int numSample, int const numChannel);

void free_hostMemory(struct hostPtr &hpp);
void free_cudaMemory(struct devicePtr &dpp);

void Fourier_beamforming_GPU_LRI(struct devicePtr &dpp, struct hostPtr hpp);

void Fourier_beamforming_GPU_HRI(struct devicePtr &dpp, struct hostPtr hpp);

void temporalFourier(Complex *d_RFIn, Complex *d_filter, int *d_coef, int numSample, int numChannel,
	cufftHandle planFFT_space, cudaStream_t stream);

void polarAxialshift(Complex *d_RFIn, int numSample, int numChannel, const float d_angle,const float d_t0, cudaStream_t stream);

void polarLateralshift(Complex *d_RFIn, float *d_kx, int numSample, int numChannel, cudaStream_t stream);

void spatialFourier(Complex *d_RFIn, Complex *d_RFIn_tranpose, int numSample, int numChannel,
	cufftHandle planFFT_space, cudaStream_t stream);

void warpFourierSpace(Complex *d_LRI, Complex *d_RFIn, cudaMemcpy3DParms myparms, const float angle, float *d_kx, int numSample, int numChannel,
	cudaStream_t stream, int streamId);

void InverseFFT(Complex *d_RFIn, Complex *d_RFIn_tranpose, int numSample, int numChannel,
	cufftHandle planFFT_time, cufftHandle planFFT_space, cudaStream_t stream);

