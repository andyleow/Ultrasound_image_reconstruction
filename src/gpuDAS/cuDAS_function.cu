#include "cuDAS.cuh"

// Declare Texture reference for 2D float texture
texture <float2, cudaTextureType1DLayered, cudaReadModeElementType> d_RF_Layertex;

__constant__ float d_t2wl;
__constant__ float d_apoConst;
__constant__ float d_ElePosC[128];
__constant__ float d_senseCutoff;

void init_hostPtr(struct hostPtr &hpp, double *h_angle, double *h_ElePos, double *h_filter, double *h_pixelMapX, double *h_pixelMapZ, double *h_t0)
{

	int i;

	hpp.ElePos = (float *)malloc(sizeof(float)*hpp.numChannel);
	hpp.angle = (float *)malloc(sizeof(float)*hpp.numAngle);
	hpp.pixelMapX = (float *)malloc(sizeof(float)*hpp.imageWidth);
	hpp.pixelMapZ = (float *)malloc(sizeof(float)*hpp.imageHeight);
	hpp.t0 = (float *)malloc(sizeof(float)*hpp.numAngle);
	hpp.filter = (float *)calloc(hpp.padSample, sizeof(float));

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
			hpp.t0[i] = ((float)h_t0[i] + fabsf(hpp.ElePos[0] * sinf(hpp.angle[i])) / hpp.c)* hpp.fs + (hpp.filterLength + 1) / 2.f;
		}
	}
	else if (hpp.focus<0.0f)
	{
		for (i = 0; i < hpp.numAngle; i++) {
			hpp.t0[i] = (float)h_t0[i] * hpp.fs + (hpp.filterLength + 1) / 2.f;  // -h_focus / h_c;  //compensate for the transmit delay
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

	hpp.init = 1;
}

//Initialise cuda memory all together
void init_cudaMemory(struct devicePtr &dpp, struct hostPtr hpp)
{
	//Select gpu
	cudaSetDevice(hpp.gpuID);

	//Initialise gpu data
	unsigned int i;

	//unsigned int sizeRF = numSample*numChannel*numAngle*numFrame * sizeof(float);
	unsigned int sizeRF = hpp.padSample*hpp.numChannel*hpp.numAngle*hpp.numFrame;
	unsigned int sizeFilt = hpp.padSample;
	unsigned int sizeIm;
    
    //Initialised cuda stream
	dpp.stream = (cudaStream_t*)malloc(nStream * sizeof(cudaStream_t));
	for (i = 0; i < nStream; i++) {
		cudaStreamCreateWithFlags(&(dpp.stream[i]), cudaStreamNonBlocking);
	}
	//cudaCheckErr();

	//Allocate pinMemory 
	cudaMallocHost((void**)&dpp.h_RFIn, sizeRF * sizeof(Complex));
	//cudaCheckErr();
	cudaMemset(dpp.h_RFIn, 0, sizeRF * sizeof(Complex));
	//cudaCheckErr();

	switch (hpp.reconMode) {
	case 0:
		sizeIm = hpp.imageHeight*hpp.imageWidth *hpp.numFrame;
		break;
	case 1:
		sizeIm = hpp.imageHeight*hpp.imageWidth *hpp.numAngle* hpp.numFrame;
		break;
	case 2:
		sizeIm = hpp.imageHeight*hpp.imageWidth *hpp.numChannel *hpp.numAngle*hpp.numFrame;
		break;
	}

	cudaMallocHost((void**)&dpp.h_LRI, sizeIm * sizeof(Complex));
	//cudaCheckErr();
	// Allocate GPU memory for each variable
	cudaMalloc((void **)& dpp.d_RFIn, sizeRF * sizeof(Complex));
	//cudaCheckErr();
	cudaMalloc((void **)& dpp.d_RFIn_temp, sizeRF * sizeof(Complex));
	//cudaCheckErr();
	cudaMalloc((void **)& dpp.d_LRI, sizeIm * sizeof(Complex));
	//cudaCheckErr();
	cudaMalloc((void **)& dpp.d_filter, sizeFilt * sizeof(Complex));
	//cudaCheckErr();
	cudaMemset(dpp.d_filter, 0, sizeFilt * sizeof(Complex));   //Filter coefficient
															   //cudaCheckErr();

	cudaMalloc((void **)& dpp.d_pixelMapX, hpp.imageWidth * sizeof(float));
	//cudaCheckErr();
	cudaMalloc((void **)& dpp.d_pixelMapZ, hpp.imageHeight * sizeof(float));
	//cudaCheckErr();

	//Copy constant variable to symbol
	float h_t2wl = hpp.fs / hpp.c;
	float h_apoConst = 3.14159265358979323f *hpp.pitch *hpp.ftx / hpp.c;  //Constant for apodization calculation (pi*d/lambda)
	cudaMemcpyToSymbol(d_apoConst, &h_apoConst, sizeof(float));
	//cudaCheckErr();
	cudaMemcpyToSymbol(d_t2wl, &h_t2wl, sizeof(float));
	//cudaCheckErr();
	cudaMemcpyToSymbol(d_ElePosC, hpp.ElePos, sizeof(float)*hpp.numChannel);
	//cudaCheckErr();
	cudaMemcpyToSymbol(d_senseCutoff, &hpp.senseCutoff, sizeof(float));
	//cudaCheckErr();

	cudaMemcpy(dpp.d_pixelMapX, hpp.pixelMapX, hpp.imageWidth * sizeof(float), cudaMemcpyHostToDevice);
	//cudaCheckErr();
	cudaMemcpy(dpp.d_pixelMapZ, hpp.pixelMapZ, hpp.imageHeight * sizeof(float), cudaMemcpyHostToDevice);
	//cudaCheckErr();
	cudaMemcpy2D(dpp.d_filter, sizeof(Complex), hpp.filter, sizeof(float), sizeof(float), hpp.padSample, cudaMemcpyHostToDevice);
	//cudaCheckErr();

	//Define hilbert coefficient
	int *h_coef;
	h_coef = (int *)calloc(hpp.padSample, sizeof(int));
	cudaMalloc((void **)&dpp.d_coef, hpp.padSample * sizeof(int));

	if (hpp.padSample % 2 == 0) {           //Even
		h_coef[0] = 1;
		h_coef[hpp.padSample / 2] = 1;
		for (i = 1; i<hpp.padSample / 2; i++)
			h_coef[i] = 2;
	}
	else {								//Odd
		h_coef[0] = 1;
		for (i = 1; i<(hpp.padSample + 1) / 2; i++)
			h_coef[i] = 2;
	}
	cudaMemcpy(dpp.d_coef, h_coef, hpp.padSample * sizeof(int), cudaMemcpyHostToDevice);
	//cudaCheckErr();

	//Initialise texture data
	cudaExtent extentDesc = make_cudaExtent(hpp.padSample, 0, nStream*hpp.numChannel);
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
	cudaMalloc3DArray(&dpp.d_RFArrayIn, &channelDesc, extentDesc, cudaArrayLayered);
	//cudaCheckErr();

	// Set texture parameters
	d_RF_Layertex.addressMode[0] = cudaAddressModeBorder;
	d_RF_Layertex.filterMode = cudaFilterModeLinear;
	d_RF_Layertex.normalized = false;    // access with normalized texture coordinates
										 // Bind the array to the texture
	cudaBindTextureToArray(d_RF_Layertex, dpp.d_RFArrayIn, d_RF_Layertex.channelDesc);
	//cudaCheckErr();

	dpp.myparms.dstArray = dpp.d_RFArrayIn;
	dpp.myparms.extent = make_cudaExtent(hpp.padSample, 1, hpp.numChannel);
	dpp.myparms.srcPos = make_cudaPos(0, 0, 0);
	dpp.myparms.kind = cudaMemcpyDeviceToDevice;

	//create batched 1D fft plan
	int batch = hpp.numChannel;
	int rank = 1;                           // --- 1D FFTs
	int n[] = { hpp.padSample };                 // --- Size of the Fourier transform
	int istride = 1, ostride = 1;           // --- Distance between two successive input/output elements
	int idist = hpp.padSample, odist = hpp.padSample; // --- Distance between batches

													  //FFT Plan for signal
	if (cufftPlanMany(&dpp.planFFT_time, rank, n,
		NULL, istride, idist,
		NULL, ostride, odist, CUFFT_C2C, batch) != CUFFT_SUCCESS) {
		printf("CUFFT error: Plan creation failed\n");
	}
	//cudaCheckErr();

	//FFT Plan for filter
	if (cufftPlanMany(&dpp.planFFT_filter, rank, n,
		NULL, istride, idist,
		NULL, ostride, odist, CUFFT_C2C, 1) != CUFFT_SUCCESS) {
		printf("CUFFT error: Plan creation failed\n");
	}
	//cudaCheckErr();

	//Execute fft for filter
	cufftSetStream(dpp.planFFT_filter, dpp.stream[0]);
	if (cufftExecC2C(dpp.planFFT_filter, (cufftComplex *)dpp.d_filter, (cufftComplex *)dpp.d_filter, CUFFT_FORWARD)) {
		printf("CUFFT Error: Unable to execute forward fft plan2\n");
	}
	cudaCheckErr();

	//cudaDeviceSynchronize();
	free(h_coef);
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
	cudaFreeHost(dpp.h_RFIn);
	cudaFreeHost(dpp.h_LRI);

	//free device memory 
	cudaFree(dpp.d_RFIn);
	cudaFree(dpp.d_RFIn_temp);
	cudaFree(dpp.d_LRI);
	cudaFree(dpp.d_filter);
	cudaFree(dpp.d_pixelMapX);
	cudaFree(dpp.d_pixelMapZ);
	cudaFree(dpp.d_coef);

	//Release memory
	cudaUnbindTexture(d_RF_Layertex);
	cudaFreeArray(dpp.d_RFArrayIn);

	//Destroy stream
	for (int i = 0; i < nStream; i++) {
		cudaStreamDestroy(dpp.stream[i]);
	}

	cufftDestroy(dpp.planFFT_time);
	cufftDestroy(dpp.planFFT_filter);
}


void firHilbert(Complex *d_RFIn, Complex *d_RFIn_temp, Complex *d_filter, int *d_coef, int numSample, int numChannel, cufftHandle planFFT_time, cudaStream_t stream)
{

	int batch = numChannel;
	cufftSetStream(planFFT_time, stream);

	//printf("Transforming signal cufftExecC2C\n");
	if (cufftExecC2C(planFFT_time, d_RFIn, d_RFIn_temp, CUFFT_FORWARD)) {
		printf("CUFFT Error: Unable to execute forward fft plan1\n"); return;
	}

	dim3 dimBlock(1024, 1, 1); //Remember max thread is only 1024
	dim3 dimGrid((numSample*batch + dimBlock.x - 1) / dimBlock.x, 1);

	//printf("Launching ComplexPointwiseMulAndScale<<< >>>\n");
	ComplexPointwiseMulAndScale << <dimGrid, dimBlock, 0, stream >> >(d_RFIn_temp, d_filter, numSample*batch, numSample, 1.0f / numSample);
	//cudaCheckErr();

	//printf("dimBlock=[%d,%d,%d]/n dimGrid=[%d,%d,%d]",dimBlock.x,dimBlock.y,dimBlock.z,dimGrid.x,dimGrid.y,dimGrid.z);
	//Peform hilbert
	EleWiseProduct << <dimGrid, dimBlock, 0, stream >> >(d_RFIn_temp, d_coef, numSample, batch);
	//cudaCheckErr();

	//Execute Reverse Transform
	if (cufftExecC2C(planFFT_time, (cufftComplex *)d_RFIn_temp, (cufftComplex *)d_RFIn, CUFFT_INVERSE) != CUFFT_SUCCESS) {
		printf("CUFFT Error: Unable to execute cufft_inverse plan\n"); return;
	}
	//cudaCheckErr();
}


__global__ void kernel_DAS_PW(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX, const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel, const float t0, const unsigned int sID)
{
	//Define linear index;
	const int ind = blockIdx.x *blockDim.x + threadIdx.x;

	if (ind< imageWidth*imageHeight)
	{
		const int j = ind / imageHeight;
		const int i = ind - imageHeight*j;
		int ChannelInd = blockIdx.y;

		float zIm = d_pixelMapZ[i];
		float xIm = d_pixelMapX[j];
		float xObj = d_pixelMapX[j] - d_ElePosC[ChannelInd];

		float t_tx = zIm * __cosf(angle) + xIm * __sinf(angle);
		float t_rx = hypotf(zIm, xObj);
		float t_bf = (t_tx + t_rx)*d_t2wl + t0;

		//Texture interpolation
		float2 temp = tex1DLayered(d_RF_Layertex, t_bf, ChannelInd + sID*numChannel);
           
		float theta = atanf(xObj/zIm - angle);
        	float  k= d_apoConst *__sinf(theta);
		float kp = __cosf(theta)/ k;
		float apo = (k == 0.0f) ? 1.0f : (__sinf(k)*kp);

		if (apo > d_senseCutoff)
		{
			d_LRI[ind].x += temp.x*apo;
			d_LRI[ind].y += temp.y*apo;
		}
	}
}


__global__ void kernel_DAS_Diverge(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX, const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel,
	const float t0, const float focus, const unsigned int sID)
{
	//Define linear index;
	const int ind = blockIdx.x *blockDim.x + threadIdx.x;

	if (ind< imageWidth*imageHeight)
	{
		
		const int j = ind / imageHeight;
		const int i = ind - imageHeight*j;
		int ChannelInd = blockIdx.y;

		//Find virtual point source
		float FocalPt_x = focus * __sinf(angle);    //coordinates of focal point changes for each angle 
		float FocalPt_z = focus * __cosf(angle);

		float zIm = d_pixelMapZ[i];
		float xIm = d_pixelMapX[j];
		float xObj = d_ElePosC[ChannelInd] - xIm;

		float t_tx = hypotf(FocalPt_x - xIm, FocalPt_z - zIm) + FocalPt_z;
		float t_rx = hypotf(zIm, xObj);
		float t_bf = (t_tx + t_rx)*d_t2wl + t0;

		//Texture interpolation
		float2 temp = tex1DLayered(d_RF_Layertex, t_bf, ChannelInd + sID*numChannel);
        
        
        	float theta = atanf(xObj/zIm - angle);
        	float  k= d_apoConst *__sinf(theta);
		float kp = __cosf(theta)/ k;
		float apo = (k == 0.0f) ? 1.0f : (__sinf(k)*kp);

		if (apo > d_senseCutoff)
		{
			d_LRI[ind].x += temp.x*apo;
			d_LRI[ind].y += temp.y*apo;
		}
	}
}


__global__ void kernel_DAS_PW_Channel(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX, const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel, const float t0, const unsigned int sID)
{
	//Define linear index;
	const int ind = blockIdx.x *blockDim.x + threadIdx.x;

	if (ind< imageWidth*imageHeight)
	{
		const int j = ind / imageHeight;
		const int i = ind - imageHeight*j;
		int ChannelInd = blockIdx.y;

		float zIm = d_pixelMapZ[i];
		float xIm = d_pixelMapX[j];
		float xObj = d_pixelMapX[j] - d_ElePosC[ChannelInd];

		float t_tx = zIm * __cosf(angle) + xIm * __sinf(angle);
		float t_rx = hypotf(zIm, xObj);
		float t_bf = (t_tx + t_rx)*d_t2wl + t0;

		//Texture interpolation
		float2 temp = tex1DLayered(d_RF_Layertex, t_bf, ChannelInd + sID*numChannel);

		float theta = atanf(xObj/zIm - angle);
        	float  k= d_apoConst *__sinf(theta);
		float kp = __cosf(theta)/ k;
		float apo = (k == 0.0f) ? 1.0f : (__sinf(k)*kp);

		if (apo > d_senseCutoff)
		{
			d_LRI[ind + ChannelInd*imageWidth*imageHeight].x += temp.x*apo;
			d_LRI[ind + ChannelInd*imageWidth*imageHeight].y += temp.y*apo;
		}
	}
}


__global__ void kernel_DAS_Diverge_Channel(float2* __restrict__  d_LRI,
	const float* __restrict__ d_pixelMapX, const float* __restrict__ d_pixelMapZ, const float angle,
	const unsigned int imageHeight, const unsigned int imageWidth, const unsigned int numChannel,
	const float t0, const float focus, const unsigned int sID)
{
	//Define linear index;
	const int ind = blockIdx.x *blockDim.x + threadIdx.x;

	if (ind< imageWidth*imageHeight)
	{
		const int j = ind / imageHeight;
		const int i = ind - imageHeight*j;
		int ChannelInd = blockIdx.y;

		//Find virtual point source
		float FocalPt_x = focus * __sinf(angle);    //coordinates of focal point changes for each angle 
		float FocalPt_z = focus * __cosf(angle);

		float zIm = d_pixelMapZ[i];
		float xIm = d_pixelMapX[j];
		float xObj = d_ElePosC[ChannelInd] - xIm;

		float t_tx = hypotf(FocalPt_x - xIm, FocalPt_z - zIm) + FocalPt_z;
		float t_rx = hypotf(zIm, xObj);
		float t_bf = (t_tx + t_rx)*d_t2wl + t0;

		//Texture interpolation
		float2 temp = tex1DLayered(d_RF_Layertex, t_bf, ChannelInd + sID*numChannel);

		float theta = atanf(xObj/zIm - angle);
        	float  k= d_apoConst *__sinf(theta);
		float kp = __cosf(theta)/ k;
		float apo = (k == 0.0f) ? 1.0f : (__sinf(k)*kp);

		if (apo > d_senseCutoff)
		{
			d_LRI[ind + ChannelInd*imageWidth*imageHeight].x += temp.x*apo;
			d_LRI[ind + ChannelInd*imageWidth*imageHeight].y += temp.y*apo;
		}
	}
}


void DAS_beamforming_GPU_LRI(struct devicePtr &dpp, struct hostPtr hpp)
{
	int i, streamId, angleId;
	int samplePerFrame = hpp.padSample*hpp.numChannel;
	int offSetIn, offSetOut;

	for (i = 0; i < hpp.numAngle*hpp.numFrame; i++)
	{
		streamId = i%nStream;
		angleId = i%hpp.numAngle;
		offSetIn = i*hpp.padSample*hpp.numChannel;
		offSetOut = i*hpp.imageHeight*hpp.imageWidth;

		//printf("Beamforming image %d", i);
		cudaStreamSynchronize(dpp.stream[streamId]);
		//Zeroing input and output memory
		cudaMemsetAsync(dpp.d_RFIn + offSetIn, 0, sizeof(Complex)*samplePerFrame, dpp.stream[streamId]);
		cudaMemsetAsync(dpp.d_LRI + offSetOut, 0, sizeof(Complex)*hpp.imageHeight*hpp.imageWidth, dpp.stream[streamId]);
		//cudaCheckErr();

		//Copy rf into device
		cudaMemcpyAsync(dpp.d_RFIn + offSetIn, dpp.h_RFIn + offSetIn, samplePerFrame * sizeof(Complex), cudaMemcpyHostToDevice, dpp.stream[streamId]);

		//1d-temporal forward FFT + fir+ Hilbert Transform
		firHilbert(dpp.d_RFIn + offSetIn, dpp.d_RFIn_temp + offSetIn, dpp.d_filter, dpp.d_coef, hpp.padSample, hpp.numChannel, dpp.planFFT_time, dpp.stream[streamId]);

		dim3 dimBlock(1024, 1, 1); //Remember max thread is only 1024
		dim3 dimGrid((hpp.imageHeight*hpp.imageWidth + dimBlock.x - 1) / dimBlock.x, hpp.numChannel, 1);

		//Copy data to texture layer
		dpp.myparms.srcPtr = make_cudaPitchedPtr(dpp.d_RFIn + offSetIn, hpp.padSample * sizeof(Complex), hpp.padSample, 1);
		dpp.myparms.dstPos = make_cudaPos(0, 0, hpp.numChannel*streamId);  //position must be between 0-2048
																		   //cudaMemcpy3D(&myparms);
		cudaMemcpy3DAsync(&dpp.myparms, dpp.stream[streamId]);

		if (hpp.focus == 0.0f)
			kernel_DAS_PW << <dimGrid, dimBlock, 0, dpp.stream[streamId] >> >(dpp.d_LRI + offSetOut, dpp.d_pixelMapX, dpp.d_pixelMapZ, hpp.angle[angleId], hpp.imageHeight, hpp.imageWidth, hpp.numChannel, hpp.t0[angleId], streamId);
		else if (hpp.focus < 0.0f)
			kernel_DAS_Diverge << <dimGrid, dimBlock, 0, dpp.stream[streamId] >> >(dpp.d_LRI + offSetOut, dpp.d_pixelMapX, dpp.d_pixelMapZ, hpp.angle[angleId], hpp.imageHeight, hpp.imageWidth, hpp.numChannel, hpp.t0[angleId], hpp.focus, streamId);

		cudaMemcpyAsync(dpp.h_LRI + offSetOut, dpp.d_LRI + offSetOut, hpp.imageHeight*hpp.imageWidth * sizeof(Complex), cudaMemcpyDeviceToHost, dpp.stream[streamId]);

	}
	cudaCheckErr();
}


void DAS_beamforming_GPU_HRI(struct devicePtr &dpp, struct hostPtr hpp)
{
	int i, streamId, angleId;
	int samplePerFrame = hpp.padSample*hpp.numChannel;
	int offSetIn, offSetOut;

	for (i = 0; i < hpp.numAngle*hpp.numFrame; i++)
	{
		streamId = i%nStream;
		angleId = i%hpp.numAngle;
		offSetIn = i*hpp.padSample*hpp.numChannel;
		offSetOut = (i / hpp.numAngle)*hpp.imageHeight*hpp.imageWidth;

		//printf("Beamforming image %d", i);

		cudaStreamSynchronize(dpp.stream[streamId]);

		//Zeroing input and output memory
		cudaMemsetAsync(dpp.d_RFIn + offSetIn, 0, sizeof(Complex)*samplePerFrame, dpp.stream[streamId]);
		if (angleId == 0)
			cudaMemsetAsync(dpp.d_LRI + offSetOut, 0, sizeof(Complex)*hpp.imageHeight*hpp.imageWidth, dpp.stream[streamId]);

		//Copy rf into device
		cudaMemcpyAsync(dpp.d_RFIn + offSetIn, dpp.h_RFIn + offSetIn, samplePerFrame * sizeof(Complex), cudaMemcpyHostToDevice, dpp.stream[streamId]);

		//1d-temporal forward FFT + fir+ Hilbert Transform
		firHilbert(dpp.d_RFIn + offSetIn, dpp.d_RFIn_temp + offSetIn, dpp.d_filter, dpp.d_coef, hpp.padSample, hpp.numChannel, dpp.planFFT_time, dpp.stream[streamId]);

		dim3 dimBlock(1024, 1, 1); //Remember max thread is only 1024
		dim3 dimGrid((hpp.imageHeight*hpp.imageWidth + dimBlock.x - 1) / dimBlock.x, hpp.numChannel, 1);

		//Copy data to texture layer
		dpp.myparms.srcPtr = make_cudaPitchedPtr(dpp.d_RFIn + offSetIn, hpp.padSample * sizeof(Complex), hpp.padSample, 1);
		dpp.myparms.dstPos = make_cudaPos(0, 0, hpp.numChannel*streamId);  //position must be between 0-2048
																		   //cudaMemcpy3D(&myparms);
		cudaMemcpy3DAsync(&dpp.myparms, dpp.stream[streamId]);

		if (hpp.focus == 0.0f)
			kernel_DAS_PW << <dimGrid, dimBlock, 0, dpp.stream[streamId] >> >(dpp.d_LRI + offSetOut, dpp.d_pixelMapX, dpp.d_pixelMapZ, hpp.angle[angleId], hpp.imageHeight, hpp.imageWidth, hpp.numChannel, hpp.t0[angleId], streamId);
		else if (hpp.focus < 0.0f)
			kernel_DAS_Diverge << <dimGrid, dimBlock, 0, dpp.stream[streamId] >> >(dpp.d_LRI + offSetOut, dpp.d_pixelMapX, dpp.d_pixelMapZ, hpp.angle[angleId], hpp.imageHeight, hpp.imageWidth, hpp.numChannel, hpp.t0[angleId], hpp.focus, streamId);

		if (angleId == (hpp.numAngle - 1)) {
			cudaMemcpyAsync(dpp.h_LRI + offSetOut, dpp.d_LRI + offSetOut, hpp.imageHeight*hpp.imageWidth * sizeof(Complex), cudaMemcpyDeviceToHost, dpp.stream[streamId]);
		}
	}
	cudaCheckErr();
}


void DAS_beamforming_GPU_CRI(struct devicePtr &dpp, struct hostPtr hpp)
{
	int i, streamId, angleId;
	int samplePerFrame = hpp.padSample*hpp.numChannel;
	int offSetIn, offSetOut;

	for (i = 0; i <hpp.numAngle*hpp.numFrame; i++)
	{
		streamId = i%nStream;
		angleId = i%hpp.numAngle;
		offSetIn = i*hpp.padSample*hpp.numChannel;
		offSetOut = i*hpp.imageHeight*hpp.imageWidth*hpp.numChannel;

		//printf("Beamforming image %d", i);
		cudaStreamSynchronize(dpp.stream[streamId]);

		//Zeroing input and output memory
		cudaMemsetAsync(dpp.d_RFIn + offSetIn, 0, sizeof(Complex)*samplePerFrame, dpp.stream[streamId]);
		cudaMemsetAsync(dpp.d_LRI + offSetOut, 0, sizeof(Complex)*hpp.imageHeight*hpp.imageWidth, dpp.stream[streamId]);
		//cudaCheckErr();

		//Copy rf into device
		cudaMemcpyAsync(dpp.d_RFIn + offSetIn, dpp.h_RFIn + offSetIn, samplePerFrame * sizeof(Complex), cudaMemcpyHostToDevice, dpp.stream[streamId]);
		//cudaCheckErr();

		//1d-temporal forward FFT + fir+ Hilbert Transform
		firHilbert(dpp.d_RFIn + offSetIn, dpp.d_RFIn_temp + offSetIn, dpp.d_filter, dpp.d_coef, hpp.padSample, hpp.numChannel, dpp.planFFT_time, dpp.stream[streamId]);

		//printf("Launching ComplexPointwiseMulAndScale<<< >>>\n");
		dim3 dimBlock(1024, 1, 1); //Remember max thread is only 1024
		dim3 dimGrid((hpp.imageHeight*hpp.imageWidth + dimBlock.x - 1) / dimBlock.x, hpp.numChannel, 1);

		//Copy data to texture layer
		dpp.myparms.srcPtr = make_cudaPitchedPtr(dpp.d_RFIn + offSetIn, hpp.padSample * sizeof(Complex), hpp.padSample, 1);
		dpp.myparms.dstPos = make_cudaPos(0, 0, hpp.numChannel*streamId);  //position must be between 0-2048
																		   //cudaMemcpy3D(&myparms);
		cudaMemcpy3DAsync(&dpp.myparms, dpp.stream[streamId]);

		if (hpp.focus == 0.0f)
			kernel_DAS_PW_Channel << <dimGrid, dimBlock, 0, dpp.stream[streamId] >> >(dpp.d_LRI + offSetOut, dpp.d_pixelMapX, dpp.d_pixelMapZ, hpp.angle[angleId], hpp.imageHeight, hpp.imageWidth, hpp.numChannel, hpp.t0[angleId], streamId);
		else if (hpp.focus < 0.0f)
			kernel_DAS_Diverge_Channel << <dimGrid, dimBlock, 0, dpp.stream[streamId] >> >(dpp.d_LRI + offSetOut, dpp.d_pixelMapX, dpp.d_pixelMapZ, hpp.angle[angleId], hpp.imageHeight, hpp.imageWidth, hpp.numChannel, hpp.t0[angleId], hpp.focus, streamId);

		cudaMemcpyAsync(dpp.h_LRI + offSetOut, dpp.d_LRI + offSetOut, hpp.imageHeight*hpp.imageWidth *hpp.numChannel * sizeof(Complex), cudaMemcpyDeviceToHost, dpp.stream[streamId]);

	}
	cudaCheckErr();
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


__global__ void EleWiseProduct(float2 *d_RFIn, int *d_coef, unsigned int numSample, unsigned int batch)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x; //row
	int ind = i%numSample;

	if (i< (numSample * batch)) {
		d_RFIn[i].x *= d_coef[ind];
		d_RFIn[i].y *= d_coef[ind];
	}

}
