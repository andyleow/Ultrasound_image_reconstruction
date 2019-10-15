#include "cuFBF.cuh"

//Define for transpose operation
#define TILE_DIM    16
#define BLOCK_ROWS  16

// Declare Texture reference for 2D float texture
texture <Complex, cudaTextureType1DLayered, cudaReadModeElementType> d_RF_Layertex;

__constant__ float d_lambda;
__constant__ float d_apoConst;
__constant__ float d_c;
__constant__ float d_pitch;
__constant__ float d_fs;
__constant__ float d_deltaX;

void init_hostPtr(struct hostPtr &hpp, double *h_angle, double *h_ElePos, double *h_filter, double *h_pixelMapX, double *h_pixelMapZ, double *h_t0)
{
	int i;

	hpp.ElePos = (float *)malloc(sizeof(float)*hpp.numChannel);
	hpp.angle = (float *)malloc(sizeof(float)*hpp.numAngle);
	hpp.pixelMapX = (float *)malloc(sizeof(float)*hpp.imageWidth);
	hpp.pixelMapZ = (float *)malloc(sizeof(float)*hpp.imageHeight);
	hpp.t0 = (float *)malloc(sizeof(float)*hpp.numAngle);
	hpp.filter = (float *)calloc(hpp.padLength, sizeof(float));

	//Preset settings
	for (i = 0; i<hpp.numChannel; i++) {
		hpp.ElePos[i] = (float)h_ElePos[i];
	}
	hpp.pitch = hpp.ElePos[1] - hpp.ElePos[0];

	for (i = 0; i<hpp.numAngle; i++) {
		hpp.angle[i] = (float)h_angle[i];
	}

	if (hpp.focus == 0.0f)
	{
		for (i = 0; i < hpp.numAngle; i++) {
			hpp.t0[i] = (float)h_t0[i] + (hpp.filterLength) / (2.f * hpp.fs);
		}
	}
	else if (hpp.focus<0.0f)
	{
		for (i = 0; i < hpp.numAngle; i++) {
			hpp.t0[i] = (float)h_t0[i] + (hpp.filterLength) / (2.f * hpp.fs); // -h_focus / h_c;  //compensate for the transmit delay
		}
	}

	for (i = 0; i<hpp.imageWidth; i++) {
		hpp.pixelMapX[i] = (float)h_pixelMapX[i];
	}

	for (i = 0; i<hpp.imageHeight; i++) {
		hpp.pixelMapZ[i] = (float)h_pixelMapZ[i];
	}

	for (i = 0; i<hpp.filterLength; i++) {
		hpp.filter[i] = (float)h_filter[i];
	}

	hpp.deltaX = hpp.pixelMapX[0] - hpp.ElePos[0];
	hpp.deltaZ = (int) roundf(hpp.pixelMapZ[0] * 2 * hpp.fs / hpp.c);
	hpp.init = 1;
}

//Initialise cuda memory all together
void init_cudaMemory(struct devicePtr &dpp, struct hostPtr hpp)
{

	//Select gpu
	cudaSetDevice(hpp.gpuID);

	int i;
	unsigned int sizeRF = hpp.padLength*hpp.padWidth*hpp.numAngle*hpp.numFrame;
	unsigned int sizeIm;
	unsigned int sizeFilt = hpp.padLength;

    //Initialised cuda stream
	dpp.stream = (cudaStream_t*)malloc(nStream * sizeof(cudaStream_t));

	for (i = 0; i < nStream; i++) {
		//cudaStreamCreate(&(stream[i]));
		cudaStreamCreateWithFlags(&(dpp.stream[i]), cudaStreamNonBlocking);
	}

	//Allocate pinned memory 
	cudaMallocHost((void**)&dpp.h_RFIn_f, sizeof(Complex)*sizeRF);
	cudaMemset(dpp.h_RFIn_f, 0, sizeof(Complex)*sizeRF);

	switch (hpp.reconMode) {
	case 0:
		sizeIm = hpp.padLength*hpp.padWidth *hpp.numFrame;
		break;
	case 1:
		sizeIm = hpp.padLength*hpp.padWidth *hpp.numAngle* hpp.numFrame;
		break;
	}
	cudaMallocHost((void**)&dpp.h_LRI_f, sizeof(Complex)*sizeIm);

	// Allocate GPU memory for each variable
	cudaMalloc((void **)& dpp.d_RFIn, sizeof(Complex)*sizeRF);
	cudaMalloc((void **)& dpp.d_RFIn_tranpose, sizeof(Complex)*sizeRF);
	cudaMalloc((void **)& dpp.d_LRI, sizeIm * sizeof(Complex));
	cudaMalloc((void **)& dpp.d_LRI_transpose, sizeIm * sizeof(Complex));

	cudaMalloc((void **)& dpp.d_filter, sizeof(Complex)*sizeFilt);
	cudaMemset(dpp.d_filter, 0, sizeof(Complex)*sizeFilt);   //Filter coefficient

	int *h_coef;
	cudaMalloc((void **)&dpp.d_coef, hpp.padLength * sizeof(int));	// Hilbert coefficient
	h_coef = (int *)calloc(hpp.padLength, sizeof(int));
	if (hpp.padLength % 2 == 0) {           //Even
		h_coef[0] = 1;
		h_coef[hpp.padLength / 2] = 1;
		for (i = 1; i<hpp.padLength / 2; i++)
			h_coef[i] = 2;
	}
	else {								//Odd
		h_coef[0] = 1;
		for (i = 1; i<(hpp.padLength + 1) / 2; i++)
			h_coef[i] = 2;
	}
	cudaMemcpy(dpp.d_coef, h_coef, hpp.padLength * sizeof(int), cudaMemcpyHostToDevice);

	
	//Predefine spatial frequency grid
	float *h_kx;
	h_kx = (float *)malloc(sizeof(float)*hpp.padWidth);
	cudaMalloc((void **)&dpp.d_kx, sizeof(float)*hpp.padWidth);
	for (i = 0; i < hpp.padWidth; i++)
	{
		if (i < hpp.padWidth / 2 + 1)
			h_kx[i] = i / hpp.pitch / hpp.padWidth;
		else
			h_kx[i] = (-hpp.padWidth + i) / hpp.pitch / hpp.padWidth;
	}

	cudaMemcpy(dpp.d_kx, h_kx, sizeof(float)*hpp.padWidth, cudaMemcpyHostToDevice);

	//Initialise texture data
	cudaExtent extentDesc = make_cudaExtent(hpp.padLength, 0, nStream*hpp.padWidth);  //Length reduce by half in analytic Fourier space
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
	cudaMalloc3DArray(&dpp.d_RFArrayIn, &channelDesc, extentDesc, cudaArrayLayered);
	//cudaMalloc3DArray(&d_RFArrayIn, &d_RF_Layertex.channelDesc, extentDesc, cudaArrayLayered);

	// Set texture parameters
	d_RF_Layertex.addressMode[0] = cudaAddressModeBorder;
	d_RF_Layertex.filterMode = cudaFilterModeLinear;
	d_RF_Layertex.normalized = false;    // access with normalized texture coordinates

										 // Bind the array to the texture
	cudaBindTextureToArray(d_RF_Layertex, dpp.d_RFArrayIn, d_RF_Layertex.channelDesc);

	dpp.myparms.dstArray = dpp.d_RFArrayIn;
	dpp.myparms.extent = make_cudaExtent(hpp.padLength, 1, hpp.padWidth);
	dpp.myparms.srcPos = make_cudaPos(0, 0, 0);
	dpp.myparms.kind = cudaMemcpyDeviceToDevice;

	//Copy constant variable to symbol
	float lambda = hpp.c / hpp.fs;
	float h_apoConst = 2 * 3.14159265358979323f *hpp.pitch;  //Constant for apodization calculation (pi*d/lambda)
	cudaMemcpyToSymbol(d_pitch, &hpp.pitch, sizeof(float));
	cudaMemcpyToSymbol(d_c, &hpp.c, sizeof(float));
	cudaMemcpyToSymbol(d_fs, &hpp.fs, sizeof(float));
	cudaMemcpyToSymbol(d_lambda, &lambda, sizeof(float));
	cudaMemcpyToSymbol(d_deltaX, &hpp.deltaX, sizeof(float));
	cudaMemcpyToSymbol(d_apoConst, &h_apoConst, sizeof(float));

	//FFT PLAN
	//Forward FFT for signal-time
	int batch = hpp.numChannel;
	int rank = 1;                           // --- 1D FFTs
	int n[] = { hpp.padLength };                 // --- Size of the Fourier transform
	int istride = 1, ostride = 1;           // --- Distance between two successive input/output elements
	int idist = hpp.padLength, odist = hpp.padLength; // --- Distance between batches

	if (cufftPlanMany(&dpp.planFFT_time_fwd, rank, n,
		NULL, istride, idist,
		NULL, ostride, odist, CUFFT_C2C, batch) != CUFFT_SUCCESS) {
		printf("CUFFT error: Plan 1 (Forward fft) creation failed \n");
		//fprintf(stderr, "CUFFT error: Plan 1 (Forward fft) creation failed \n"); return;
	}

	//Reverse FFT for signal-time
	batch = hpp.padWidth;
	n[0] = hpp.padLength;                 // --- Size of the Fourier transform
	idist = hpp.padLength;
	odist = hpp.padLength; // --- Distance between batches

	if (cufftPlanMany(&dpp.planFFT_time_rvs, rank, n,
		NULL, istride, idist,
		NULL, ostride, odist, CUFFT_C2C, batch) != CUFFT_SUCCESS) {
		printf("CUFFT error: Plan 1 (Forward fft) creation failed \n");
		//fprintf(stderr, "CUFFT error: Plan 1 (Forward fft) creation failed \n"); return;
	}

	//forward fft lateral space signal
	batch = hpp.padLength;
	n[0] = hpp.padWidth;                 // --- Size of the Fourier transform
	idist = hpp.padWidth;
	odist = hpp.padWidth; // --- Distance between batches
						  //FFT for signal
	if (cufftPlanMany(&dpp.planFFT_space_fwd, rank, n,
		NULL, istride, idist,
		NULL, ostride, odist, CUFFT_C2C, batch) != CUFFT_SUCCESS) {
		printf("CUFFT error: Plan creation failed\n");
	}

	//Reverse fft lateral space signal
	batch = hpp.padLength;
	n[0] = hpp.padWidth;                 // --- Size of the Fourier transform
	idist = hpp.padWidth;
	odist = hpp.padWidth; // --- Distance between batches
						  //FFT for signal
	if (cufftPlanMany(&dpp.planFFT_space_rvs, rank, n,
		NULL, istride, idist,
		NULL, ostride, odist, CUFFT_C2C, batch) != CUFFT_SUCCESS) {
		printf("CUFFT error: Plan creation failed\n");
	}

	//forward FFT for filter
	n[0] = hpp.padLength;                 // --- Size of the Fourier transform
	idist = hpp.padLength;
	odist = hpp.padLength; // --- Distance between batches

	if (cufftPlanMany(&dpp.planFFT_filter, rank, n,
		NULL, istride, idist,
		NULL, ostride, odist, CUFFT_C2C, 1) != CUFFT_SUCCESS) {
		printf("CUFFT error: Plan 2 (Forward fft) creation failed\n");
	}

	//Copy variable from host to device
	cudaMemcpy2D(dpp.d_filter, sizeof(Complex), hpp.filter, sizeof(float), sizeof(float), hpp.padLength, cudaMemcpyHostToDevice);

	if (cufftExecC2C(dpp.planFFT_filter, (cufftComplex *)dpp.d_filter, (cufftComplex *)dpp.d_filter, CUFFT_FORWARD)) {
		printf("CUFFT Error: Unable to execute forward fft plan2\n");
	}

	cudaCheckErr();
	free(h_coef);
	free(h_kx);

}

void free_hostMemory(struct hostPtr &hpp)
{
	free(hpp.ElePos);
	free(hpp.angle);
	free(hpp.pixelMapX);
	free(hpp.pixelMapZ);
	free(hpp.filter);
	hpp.init = 0;

}

void free_cudaMemory(struct devicePtr &dpp)
{
	//Free pinned hostr memory 
	cudaFreeHost(dpp.h_RFIn_f);
	cudaFreeHost(dpp.h_LRI_f);

	//free device memory 
	cudaFree(dpp.d_filter);
	cudaFree(dpp.d_RFIn);
	cudaFree(dpp.d_RFIn_tranpose);
	cudaFree(dpp.d_LRI);
	cudaFree(dpp.d_LRI_transpose);
	cudaFree(dpp.d_coef);
	cudaFree(dpp.d_kx);

	//Release memory
	cudaUnbindTexture(d_RF_Layertex);
	cudaFreeArray(dpp.d_RFArrayIn);

	//Destroy stream
	for (int i = 0; i < nStream; i++) {
		cudaStreamDestroy(dpp.stream[i]);
	}

	cufftDestroy(dpp.planFFT_filter);
	cufftDestroy(dpp.planFFT_time_fwd);
	cufftDestroy(dpp.planFFT_time_rvs);
	cufftDestroy(dpp.planFFT_space_fwd);
	cufftDestroy(dpp.planFFT_space_rvs);
}


void Fourier_beamforming_GPU_LRI(struct devicePtr &dpp, struct hostPtr hpp)
{
	int i, streamId, angleId;
	int samplePerFrame = hpp.padLength*hpp.padWidth;
	int offsetIn, offsetOut;

	for (i = 0; i<hpp.numAngle*hpp.numFrame; i++)
	{
		streamId = i%nStream;
		angleId = i%hpp.numAngle;
		offsetIn = i*samplePerFrame;
		offsetOut = i*samplePerFrame;

		cudaStreamSynchronize(dpp.stream[streamId]);

		//Copy data to device
		cudaMemsetAsync(dpp.d_RFIn + offsetIn + samplePerFrame / 2, 0, sizeof(Complex)*samplePerFrame / 2, dpp.stream[streamId]);  //Zero out half memory not init 
		cudaMemcpyAsync(dpp.d_RFIn + offsetIn, dpp.h_RFIn_f + offsetIn, samplePerFrame / 2 * sizeof(Complex), cudaMemcpyHostToDevice, dpp.stream[streamId]);
		cudaMemsetAsync(dpp.d_LRI + offsetOut, 0, sizeof(Complex)*samplePerFrame, dpp.stream[streamId]);

		//1d-temporal forward FFT + fir+ Hilbert Transform
		//printf("\nTemporal fft");
		temporalFourier(dpp.d_RFIn + offsetIn, dpp.d_filter, dpp.d_coef, hpp.padLength, hpp.numChannel, dpp.planFFT_time_fwd, dpp.stream[streamId]);


		//Polar shifting to time zero
		//printf("\npolar shifting");
		polarAxialshift(dpp.d_RFIn + offsetIn, hpp.padLength, hpp.numChannel, hpp.angle[angleId], hpp.t0[angleId], dpp.stream[streamId]);

		//Spatial Fourier Transform
		//printf("\nSpatial fft");
		spatialFourier(dpp.d_RFIn + offsetIn, dpp.d_RFIn_tranpose + offsetIn, hpp.padLength, hpp.padWidth, dpp.planFFT_space_fwd, dpp.stream[streamId]);

        //Lateral shift for extra FOV
		//printf("\Extending lateral FOV");
		polarLateralshift(dpp.d_RFIn + offsetIn, dpp.d_kx, hpp.padLength, hpp.padWidth, dpp.stream[streamId]);
		
		//Warp fourier space
		//printf("\nFourier Interpolation");
		warpFourierSpace(dpp.d_LRI + offsetOut, dpp.d_RFIn + offsetIn, dpp.myparms, hpp.angle[angleId], dpp.d_kx, hpp.padLength, hpp.padWidth, dpp.stream[streamId], streamId);

		// 2d inverse Fourier Transform
		//printf("\nInverse FFT");
		InverseFFT(dpp.d_LRI + +offsetOut, dpp.d_RFIn_tranpose + offsetIn, hpp.padLength, hpp.padWidth, dpp.planFFT_time_rvs, dpp.planFFT_space_rvs, dpp.stream[streamId]);
		
		//Copy all data back to host
		//printf("\nGPU-CPU Transfer");
		cudaMemcpyAsync(dpp.h_LRI_f + offsetIn, dpp.d_LRI + offsetOut, samplePerFrame * sizeof(Complex), cudaMemcpyDeviceToHost, dpp.stream[streamId]);
		//cudaStreamQuery(dpp.stream[streamId]);
	}
 cudaCheckErr();
}


void Fourier_beamforming_GPU_HRI(struct devicePtr &dpp, struct hostPtr hpp)
{
	int i, streamId, angleId;
	int samplePerFrame = hpp.padLength*hpp.padWidth;
	int offsetIn, offsetOut;

	for (i = 0; i<hpp.numAngle*hpp.numFrame; i++)
	{
		streamId = i%nStream;
		angleId = i%hpp.numAngle;
		offsetIn = i*samplePerFrame;
		offsetOut = (i / hpp.numAngle)*samplePerFrame;

		cudaStreamSynchronize(dpp.stream[streamId]);

		//Copy data to device
		cudaMemsetAsync(dpp.d_RFIn + offsetIn + samplePerFrame / 2, 0, sizeof(Complex)*samplePerFrame / 2, dpp.stream[streamId]);
		cudaMemcpyAsync(dpp.d_RFIn + offsetIn, dpp.h_RFIn_f + offsetIn, samplePerFrame / 2 * sizeof(Complex), cudaMemcpyHostToDevice, dpp.stream[streamId]);
		if (angleId == 0)
			cudaMemsetAsync(dpp.d_LRI + offsetOut, 0, sizeof(Complex)*samplePerFrame, dpp.stream[streamId]);


		//1d-temporal forward FFT + fir+ Hilbert Transform
		//printf("\nTemporal fft");
		temporalFourier(dpp.d_RFIn + offsetIn, dpp.d_filter, dpp.d_coef, hpp.padLength, hpp.numChannel, dpp.planFFT_time_fwd, dpp.stream[streamId]);


		//Polar shifting to time zero
		//printf("\npolar shifting");
		polarAxialshift(dpp.d_RFIn + offsetIn, hpp.padLength, hpp.numChannel, hpp.angle[angleId], hpp.t0[angleId], dpp.stream[streamId]);

		//Spatial Fourier Transform
		//printf("\nSpatial fft");
		spatialFourier(dpp.d_RFIn + offsetIn, dpp.d_RFIn_tranpose + offsetIn, hpp.padLength, hpp.padWidth, dpp.planFFT_space_fwd, dpp.stream[streamId]);

		//Lateral shift for extra FOV
		//printf("\Extending lateral FOV");
		polarLateralshift(dpp.d_RFIn + offsetIn, dpp.d_kx, hpp.padLength, hpp.padWidth, dpp.stream[streamId]);


		//Warp fourier space
		//printf("\nFourier Interpolation");
		warpFourierSpace(dpp.d_LRI + offsetOut, dpp.d_RFIn + offsetIn, dpp.myparms, hpp.angle[angleId], dpp.d_kx, hpp.padLength, hpp.padWidth, dpp.stream[streamId], streamId);

		//Copy all data back to host
		//printf("\nGPU-CPU Transfer");
		if (angleId == (hpp.numAngle - 1))
		{
			InverseFFT(dpp.d_LRI + offsetOut, dpp.d_LRI_transpose + offsetOut, hpp.padLength, hpp.padWidth, dpp.planFFT_time_rvs, dpp.planFFT_space_rvs, dpp.stream[streamId]);
			cudaMemcpyAsync(dpp.h_LRI_f + offsetOut, dpp.d_LRI + offsetOut, samplePerFrame * sizeof(Complex), cudaMemcpyDeviceToHost, dpp.stream[streamId]);
		}
	}
 cudaCheckErr();    
}



void temporalFourier(Complex *d_RFIn, Complex *d_filter, int *d_coef, int numSample, int numChannel,
	cufftHandle planFFT_time, cudaStream_t stream)
{
	//Batch size is limited due to memory contrain
	//Max batch size approximate (1-8)GB/(2*sizeof(Complex)*n[0]....*n[rank])
	int batch = numChannel;
	cufftSetStream(planFFT_time, stream);

	//printf("Transforming signal cufftExecC2C\n");
	if (cufftExecC2C(planFFT_time, d_RFIn, d_RFIn, CUFFT_FORWARD)) {
		printf("CUFFT Error: Unable to execute forward fft plan1\n"); return;
	}

	//printf("Launching ComplexPointwiseMulAndScale<<< >>>\n");
	dim3 dimBlock(1024, 1, 1); //Remember max thread is only 1024
	dim3 dimGrid((numSample*batch + dimBlock.x - 1) / dimBlock.x, 1);

	ComplexPointwiseMulAndScale << <dimGrid, dimBlock, 0, stream >> >(d_RFIn, d_filter, numSample*batch, numSample, 1.0f / (numSample*numChannel/2));

	//Peform hilbert
	//printf("dimBlock=[%d,%d,%d]/n dimGrid=[%d,%d,%d]",dimBlock.x,dimBlock.y,dimBlock.z,dimGrid.x,dimGrid.y,dimGrid.z);
	EleWiseProduct << <dimGrid, dimBlock, 0, stream >> >(d_RFIn, d_coef, numSample, batch);
}


void polarAxialshift(Complex *d_RFIn, int numSample, int numChannel, const float h_angle, const float h_t0, cudaStream_t stream)
{
	dim3 dimBlock(1024, 1, 1);
	dim3 dimGrid(((numSample / 2 + 1) + dimBlock.x - 1) / dimBlock.x, numChannel, 1);
	//printf("dimBlock=[%d,%d,%d]/n dimGrid=[%d,%d,%d]", dimBlock.x, dimBlock.y, dimBlock.z, dimGrid.x, dimGrid.y, dimGrid.z);

	kernel_complexAxialShift << <dimGrid, dimBlock, 0, stream >> >(d_RFIn, h_angle, h_t0, numSample, numChannel);
}

void polarLateralshift(Complex *d_RFIn, float *d_kx, int numSample, int numChannel, cudaStream_t stream)
{
	dim3 dimBlock(1024, 1, 1);
	dim3 dimGrid(((numSample + 1) + dimBlock.x - 1) / dimBlock.x, numChannel, 1);
	//printf("dimBlock=[%d,%d,%d]/n dimGrid=[%d,%d,%d]", dimBlock.x, dimBlock.y, dimBlock.z, dimGrid.x, dimGrid.y, dimGrid.z);

	kernel_complexLateralShift << <dimGrid, dimBlock, 0, stream >> >(d_RFIn, d_kx, numSample, numChannel);
}

void spatialFourier(Complex *d_RFIn, Complex *d_RFIn_tranpose, int numSample, int numChannel, cufftHandle planFFT_space, cudaStream_t stream)
{
	//Preforming matrix transpose, 1d fft, and matrix transpose for memory coalescent 

	//  transpose operation
	dim3 dimBlock(TILE_DIM, TILE_DIM);
	dim3 dimGrid(numSample / TILE_DIM, numChannel / TILE_DIM);


	kernel_transpose << <dimGrid, dimBlock, 0, stream >> > (d_RFIn_tranpose, d_RFIn, numSample, numChannel);

	// 1d fft 
	//cudaStreamSynchronize(stream);
	//printf("\nFFT processing frame : %d", i);
	cufftSetStream(planFFT_space, stream);

	//printf("Transforming signal cufftExecC2C\n");
	if (cufftExecC2C(planFFT_space, d_RFIn_tranpose, d_RFIn_tranpose, CUFFT_FORWARD)) {
		printf("CUFFT Error: Unable to execute forward fft plan1\n"); return;
	}

	// matrix tranpose
	dimGrid.y = numSample / TILE_DIM;
	dimGrid.x = numChannel / TILE_DIM;

	//cudaStreamSynchronize(stream);
	kernel_transpose << <dimGrid, dimBlock, 0, stream >> > (d_RFIn, d_RFIn_tranpose, numChannel, numSample);

}


void warpFourierSpace(Complex *d_LRI, Complex *d_RFIn, cudaMemcpy3DParms myparms, const float angle, float *d_kx, int numSample, int numChannel,
	cudaStream_t stream, int streamId)
{
	dim3 dimBlock(1024, 1, 1); //Remember max thread is only 1024
							   // 	dim3 dimGrid(( imageHeight  + dimBlock.x - 1) / dimBlock.x,
							   //                  (imageWidth + dimBlock.y - 1) / dimBlock.y);
	dim3 dimGrid(((numSample / 2 + 1) + dimBlock.x - 1) / dimBlock.x,
		numChannel, 1);
	//printf("dimBlock=[%d,%d,%d]/n dimGrid=[%d,%d,%d]", dimBlock.x, dimBlock.y, dimBlock.z, dimGrid.x, dimGrid.y, dimGrid.z);

	//cudaStreamSynchronize(stream);
	//Load data to texture memory
	myparms.srcPtr = make_cudaPitchedPtr(d_RFIn, numSample * sizeof(Complex), numSample, 1);
	myparms.dstPos = make_cudaPos(0, 0, numChannel*streamId);  //position must be between 0-2048
															   //cudaMemcpy3D(&myparms);
	cudaMemcpy3DAsync(&myparms, stream);
	kernel_FourierInterpolation << <dimGrid, dimBlock, 0, stream >> >(d_LRI, d_kx, angle, numSample, numChannel, streamId);
}

void InverseFFT(Complex *d_RFIn, Complex *d_RFIn_tranpose, int numSample, int numChannel, cufftHandle planFFT_time,
	cufftHandle planFFT_space, cudaStream_t stream)
{
	//Perform two 1d fft, transpose the data before doing fft on x direction 

	// 1d fft across depth direction
	//cudaStreamSynchronize(stream);
	//printf("\nFFT processing frame : %d", i);	
	cufftSetStream(planFFT_time, stream);

	//printf("Transforming signal cufftExecC2C\n");
	if (cufftExecC2C(planFFT_time, d_RFIn, d_RFIn, CUFFT_INVERSE)) {
		printf("CUFFT Error: Unable to execute forward fft plan1\n"); return;
	}

	dim3 dimBlock(TILE_DIM, TILE_DIM);
	dim3 dimGrid(numSample / TILE_DIM, numChannel / TILE_DIM);

	// transpose operation
	kernel_transpose << <dimGrid, dimBlock, 0, stream >> > (d_RFIn_tranpose, d_RFIn, numSample, numChannel);

	// 1d fft 
	//cudaStreamSynchronize(stream);
	//printf("\nFFT processing frame : %d", i);

	cufftSetStream(planFFT_space, stream);

	//printf("Transforming signal cufftExecC2C\n");
	if (cufftExecC2C(planFFT_space, d_RFIn_tranpose, d_RFIn_tranpose, CUFFT_INVERSE)) {
		printf("CUFFT Error: Unable to execute forward fft plan1\n"); return;
	}

	// matrix tranpose
	dimGrid.y = numSample / TILE_DIM;
	dimGrid.x = numChannel / TILE_DIM;

	kernel_transpose << <dimGrid, dimBlock, 0, stream >> > (d_RFIn, d_RFIn_tranpose, numChannel, numSample);
}

__global__ void kernel_transpose(Complex *odata, Complex *idata, int width, int height)
{
	__shared__ Complex tile[TILE_DIM][TILE_DIM + 1];

	unsigned int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

	if ((xIndex < width) && (yIndex < height))
	{
		unsigned int index_in = xIndex + (yIndex)*width;
		tile[threadIdx.y][threadIdx.x].x = idata[index_in].x;
		tile[threadIdx.y][threadIdx.x].y = idata[index_in].y;
	}
	__syncthreads();

	xIndex = blockIdx.y * TILE_DIM + threadIdx.x;
	yIndex = blockIdx.x * TILE_DIM + threadIdx.y;

	if ((xIndex < height) && (yIndex < width))
	{
		unsigned int index_out = xIndex + (yIndex)*height;
		odata[index_out].x = tile[threadIdx.x][threadIdx.y].x;
		odata[index_out].y = tile[threadIdx.x][threadIdx.y].y;
	}

}

static __device__ __host__ inline Complex ComplexAdd(Complex a, Complex b)
{
	Complex c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	return c;
}

static __device__ __host__ inline Complex ComplexScale(const Complex a, const float s)
{
	Complex c;
	c.x = s * a.x;
	c.y = s * a.y;
	return c;
}

// Complex multiplication
static __device__ __host__ inline Complex ComplexMul(const Complex a, const Complex b)
{
	Complex c;
	c.x = a.x * b.x - a.y * b.y;
	c.y = a.x * b.y + a.y * b.x;
	return c;
}

// Complex pointwise multiplication
static __global__ void ComplexPointwiseMulAndScale(Complex *a, const Complex *b, int size, int numSample, float scale)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	const int n = i%numSample;

	if (i < size)
	{
		a[i] = ComplexScale(ComplexMul(a[i], b[n]), scale);
	}
}


__global__ void EleWiseProduct(Complex *d_RFIn, int *d_coef, unsigned int numSample, unsigned int batch)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x; //row
	int ind = i%numSample;

	if (i< (numSample * batch)) {
		d_RFIn[i].x *= d_coef[ind];
		d_RFIn[i].y *= d_coef[ind];
	}

}

__global__ void kernel_complexAxialShift(Complex *d_RFIn, float angle, float t0, unsigned int numSample, int numChannel)
{
	//Define linear and 2d index;
	const int i = blockIdx.x *blockDim.x + threadIdx.x;
	const int j = blockIdx.y;
	const int ind = j*numSample + i;


	if (i < (numSample / 2 + 1))
	{
		const float sinA = __sinf(angle);
		const float FRAC_TWOPI_F = 6.283185307179586f*d_fs / numSample;
		const float PITCH_OVER_C = d_pitch / d_c;

		//tau=param.pitch / param.c*(sinA*(0:nx - 1) + (nx - 1)*abs(sinA)*(sinA<0)) + param.t0(i);
		float tmin = (sinA < 0) ? j - (numChannel - 1) : j;
		float t = tmin*sinA*PITCH_OVER_C + t0;
		float phi = (FRAC_TWOPI_F*i)*t;

		//Complex exponential multiplication (S =S *exp(j*2*pi*f*t))
		Complex b = { __cosf(phi),__sinf(phi) };
		d_RFIn[ind] = ComplexMul(d_RFIn[ind], b);
	}
}

__global__ void kernel_complexLateralShift(Complex *d_RFIn, float *d_kx, unsigned int numSample, int numChannel)
{
	//Define linear and 2d index;
	const int i = blockIdx.x *blockDim.x + threadIdx.x;
	const int j = blockIdx.y;
	const int ind = j*numSample + i;


	if (i < (numSample / 2 + 1))
	{
		const float TWOPI_SHIFT = 6.283185307179586f *d_deltaX;
		float phi = TWOPI_SHIFT * d_kx[j];

		//Complex exponential multiplication (S =S *exp(j*2*pi*f*t))
		Complex b = { __cosf(phi),__sinf(phi) };
		d_RFIn[ind] = ComplexMul(d_RFIn[ind], b);
	}
}

__global__ void kernel_FourierInterpolation(Complex* __restrict__  d_LRI,
	float* __restrict__ kx, const float angle, const unsigned int numSample, const unsigned int numChannel, const unsigned int sID)
{

	const int i = blockIdx.x *blockDim.x + threadIdx.x;
	const int j = blockIdx.y;


	if (i< (numSample / 2 + 1))
	{
		const int ind = j*numSample + i;
		const float nLambda = numSample *d_lambda;

		const float eps = 0.00000011920929f;
		float kzi = 2.0f* i / nLambda;
		float kxj = kx[j];
		float cosA = __cosf(angle);
		float sinA = __sinf(angle);

		//Evanescent wave weight
		float w1 = kxj * cosA / (sinA + 1);
		float w2 = kxj * cosA / (sinA - 1);
		float wEva = (kzi >= w1) * (kzi >= w2);

		//Element sensitivity
		float w3 = d_apoConst *kxj;
		float wEleSense = (w3 == 0.0f) ? 1.0f : (__sinf(w3) / w3);

		float kd = (kzi * cosA + kxj *sinA + eps);
		float oneOverkd = 1 / kd;
		float kn = kzi*kzi + kxj*kxj;

		float f = 0.5f *kn * oneOverkd *nLambda;  //Change to lambda

												  //Texture interpolation
		Complex temp = tex1DLayered(d_RF_Layertex, f, j + sID*numChannel);

		d_LRI[ind].x += temp.x *wEva *wEleSense;
		d_LRI[ind].y += temp.y *wEva *wEleSense;
	}
}
