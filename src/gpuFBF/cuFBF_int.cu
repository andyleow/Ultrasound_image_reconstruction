#include "cuFBF.cuh"

#define MODE 0
#define INPUTRF 1
#define ANGLE 2
#define ELEPOS 3
#define FILTER 4
#define PIXELMAP_X 5
#define PIXELMAP_Z 6
#define T0 7
#define FS 8
#define FTX 9
#define C 10
#define FOCUS 11
#define RECONMODE 12
#define GPUID 13

#define OUTPUTIMAGE 0

struct hostPtr hpp;
struct devicePtr dpp;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int h, i, j, nDim;
	int mode = (int)((mxGetPr(prhs[MODE]))[0]);

	int16_t *h_dataIn_r;
	double *h_angle, *h_ElePos, *h_t0;
	double *h_pixelMapX, *h_pixelMapZ;
	double *h_filter;


	switch (mode) {
	case 0:
		//Check status and execute function if memory allocated
		if (hpp.init == 0)
			mexErrMsgTxt("memory not allocated");
		else {
			if (mexIsLocked()) {
				//printf("\nmemory allocated and locked");
				break;
			}
		}

	case 1:
		//Initialise input
		if (hpp.init == 1) {
			//clear memory first
			free_hostMemory(hpp);
			free_cudaMemory(dpp);
		}

		//Read INPUT data
		h_angle = mxGetPr(prhs[ANGLE]);
		h_ElePos = mxGetPr(prhs[ELEPOS]);
		h_filter = mxGetPr(prhs[FILTER]);
		h_pixelMapX = mxGetPr(prhs[PIXELMAP_X]);
		h_pixelMapZ = mxGetPr(prhs[PIXELMAP_Z]);
		h_t0 = mxGetPr(prhs[T0]);

		nDim = mxGetNumberOfDimensions(prhs[INPUTRF]);
		hpp.numSample = ((mxGetDimensions(prhs[INPUTRF]))[0]);
		hpp.numChannel = ((mxGetDimensions(prhs[INPUTRF]))[1]);
		hpp.numAngle = 1;
		hpp.numFrame = 1;
		switch (nDim) {
		case 3:
			hpp.numAngle = ((mxGetDimensions(prhs[ANGLE]))[0] * (mxGetDimensions(prhs[ANGLE]))[1]);
			hpp.numFrame = 1;
			break;
		case 4:
            hpp.numAngle = ((mxGetDimensions(prhs[ANGLE]))[0] * (mxGetDimensions(prhs[ANGLE]))[1]);
			hpp.numFrame = ((mxGetDimensions(prhs[INPUTRF]))[3]);
		}
		//          mexPrintf("%d, %d, %d, %d", hpp.numSample,hpp.numChannel, hpp.numAngle, hpp.numFrame);

		hpp.filterLength = ((mxGetDimensions(prhs[FILTER]))[0] * (mxGetDimensions(prhs[FILTER]))[1]);

		hpp.fs = (float)((mxGetPr(prhs[FS]))[0]);
		hpp.ftx = (float)((mxGetPr(prhs[FTX]))[0]);
		hpp.c = (float)((mxGetPr(prhs[C]))[0]);
		hpp.focus = (float)((mxGetPr(prhs[FOCUS]))[0]);
		hpp.reconMode = (int)((mxGetPr(prhs[RECONMODE]))[0]);
		hpp.gpuID = (int)((mxGetPr(prhs[GPUID]))[0]);

		hpp.imageWidth = ((mxGetDimensions(prhs[PIXELMAP_X]))[0] * (mxGetDimensions(prhs[PIXELMAP_X]))[1]);
		hpp.imageHeight = ((mxGetDimensions(prhs[PIXELMAP_Z]))[0] * (mxGetDimensions(prhs[PIXELMAP_Z]))[1]);

		hpp.padLength = powf(2.0f, ceilf(log2f(hpp.numSample + 1)));
		if (hpp.imageWidth>hpp.numChannel)
			hpp.padWidth = powf(2.0f, ceilf(log2f(hpp.imageWidth + 1)));
		else
			hpp.padWidth = powf(2.0f, ceilf(log2f(hpp.numChannel + 1)));
		//         mexPrintf("\n%d, %d  %d", hpp.imageWidth, hpp.imageHeight, hpp.init);

		//Iinitialise host memory
		//          mexPrintf("Init hostPtr");
		init_hostPtr(hpp, h_angle, h_ElePos, h_filter, h_pixelMapX, h_pixelMapZ, h_t0);

		//Initialise cuda memory
		init_cudaMemory(dpp, hpp);
		mexLock();
		mexPrintf("Memory initialzed");
		
		return;
	case -1:
		mexPrintf("Mode 2");
		if (hpp.init == 1) {
			
			free_hostMemory(hpp);
			free_cudaMemory(dpp);
			mexUnlock();
			mexPrintf("Memory cleared");
		}
		else
			mexErrMsgTxt("memory not allocated");
		return;
	}


	//     mexPrintf("\nCopy data %d, %d %d %d %d", hpp.padSample, hpp.numSample, hpp.numChannel,hpp.numFrame,hpp.numAngle);
	if (mxIsInt16(prhs[INPUTRF]))
		h_dataIn_r = (int16_t *)mxGetPr(prhs[INPUTRF]);
	else
		mexErrMsgTxt("RF input datatype must be int16");

	for (h = 0; h < hpp.numAngle*hpp.numFrame; h++) {
		for (i = 0; i < hpp.numChannel; i++) {
			for (j = 0; j < hpp.numSample; j++) {
				dpp.h_RFIn_f[h*(hpp.padLength*hpp.padWidth) + i*hpp.padLength + j].x = (float)h_dataIn_r[h*(hpp.numSample*hpp.numChannel) + i*hpp.numSample + j];
			}
		}
	}

	int nFrames;
	mwSize dim[4] = { hpp.imageHeight,hpp.imageWidth,hpp.numAngle,hpp.numFrame };

	switch (hpp.reconMode) {
	case 0:

		//Dyn Apo-HRI-Images
		Fourier_beamforming_GPU_HRI(dpp, hpp);
		nFrames = hpp.numFrame;
		dim[2] = hpp.numFrame; dim[3] = 1;
		break;

	case 1:
		//Dynamic Apo-LRI-Images
		Fourier_beamforming_GPU_LRI(dpp, hpp);
		nFrames = hpp.numAngle*hpp.numFrame;
		break;
	}

	// Return output images
	plhs[OUTPUTIMAGE] = mxCreateNumericArray(4, dim, mxSINGLE_CLASS, mxCOMPLEX); //create nd single
	hpp.LRI_r = (float *)mxGetData(plhs[OUTPUTIMAGE]);
	hpp.LRI_i = (float *)mxGetImagData(plhs[OUTPUTIMAGE]);

	//Tranfer to matlab workspace
	for (h = 0; h < nFrames; h++) {
		for (i = 0; i < hpp.imageWidth; i++) {
			for (j = 0; j < hpp.imageHeight; j++) {
				hpp.LRI_r[h*hpp.imageHeight*hpp.imageWidth + i*hpp.imageHeight + j] = dpp.h_LRI_f[(h*hpp.padLength*hpp.padWidth + i*hpp.padLength +j +hpp.deltaZ)].x;
				hpp.LRI_i[h*hpp.imageHeight*hpp.imageWidth + i*hpp.imageHeight + j] = dpp.h_LRI_f[(h*hpp.padLength*hpp.padWidth + i*hpp.padLength +j +hpp.deltaZ)].y;
			}

		}
	}
}

