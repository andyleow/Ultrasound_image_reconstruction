#include "cuDAS.cuh"

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
#define ELESENSE 12
#define RECONMODE 13
#define GPUID 14

#define OUTPUTIMAGE 0

struct devicePtr dpp;
struct hostPtr hpp;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j, nDim;
	int mode = (int)((mxGetPr(prhs[MODE]))[0]);

	float *h_dataIn_r;
	double *h_angle, *h_ElePos, *h_t0;
	double *h_pixelMapX, *h_pixelMapZ;
	double *h_filter;

	switch (mode) {
	case 0:
		//Check status and execute function if memory allocated
		if (hpp.init == 0)
			mexErrMsgTxt("memory not allocated");
		else {
			if (mexIsLocked) {
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
        hpp.numAngle =1;
        hpp.numFrame = 1;
		switch (nDim) {			
		case 3:
			hpp.numAngle = ((mxGetDimensions(prhs[ANGLE]))[0] *(mxGetDimensions(prhs[ANGLE]))[1]);
			hpp.numFrame = 1;
            break;
		case 4:
			hpp.numAngle = ((mxGetDimensions(prhs[ANGLE]))[0] *(mxGetDimensions(prhs[ANGLE]))[1]);
			hpp.numFrame = ((mxGetDimensions(prhs[INPUTRF]))[3]);
		}
//         mexPrintf("%d, %d, %d, %d", hpp.numSample,hpp.numChannel, hpp.numAngle, hpp.numFrame);
      
		hpp.filterLength = ((mxGetDimensions(prhs[FILTER]))[0] * (mxGetDimensions(prhs[FILTER]))[1]);
		hpp.padSample = (int) powf(2.0f, ceilf(log2f(hpp.numSample)));

		hpp.fs = (float)((mxGetPr(prhs[FS]))[0]);
		hpp.ftx = (float)((mxGetPr(prhs[FTX]))[0]);
		hpp.c = (float)((mxGetPr(prhs[C]))[0]);
		hpp.focus = (float)((mxGetPr(prhs[FOCUS]))[0]);
		hpp.senseCutoff= (float)((mxGetPr(prhs[ELESENSE]))[0]);
		hpp.reconMode = (int)((mxGetPr(prhs[RECONMODE]))[0]);
		hpp.gpuID = (int)((mxGetPr(prhs[GPUID]))[0]);

		hpp.imageWidth = ((mxGetDimensions(prhs[PIXELMAP_X]))[0] * (mxGetDimensions(prhs[PIXELMAP_X]))[1]);
		hpp.imageHeight = ((mxGetDimensions(prhs[PIXELMAP_Z]))[0] * (mxGetDimensions(prhs[PIXELMAP_Z]))[1]);
//        mexPrintf("\n%d, %d  %d", hpp.imageWidth, hpp.imageHeight, hpp.init);

		//Iinitialise host memory
//          mexPrintf("Init hostPtr");
        init_hostPtr(hpp, h_angle, h_ElePos, h_filter, h_pixelMapX, h_pixelMapZ, h_t0);
		
		//Initialise cuda memory
//          mexPrintf("Init devicePtr");
        init_cudaMemory(dpp,hpp);
		mexLock();
//         mexPrintf("\nmex Lock");
		return;
	case -1:
// 		mexPrintf("Mode 2");
		if (hpp.init == 1) {
			//mexPrintf("clear input");
			free_hostMemory(hpp);
			free_cudaMemory(dpp);
			mexUnlock();
// 			mexPrintf("cleared");
		}
		else
			mexErrMsgTxt("memory not allocated");
		return;
	}

//     mexPrintf("\nCopy data %d, %d %d %d %d", hpp.padSample, hpp.numSample, hpp.numChannel,hpp.numFrame,hpp.numAngle);
	h_dataIn_r = (float *)mxGetPr(prhs[INPUTRF]);
	for (i = 0; i < hpp.numAngle*hpp.numFrame*hpp.numChannel; i++) {
		for (j = 0; j < hpp.numSample; j++) {
			dpp.h_RFIn[i*hpp.padSample + j].x = (float)h_dataIn_r[i*hpp.numSample + j];
		}
	}
//      mexPrintf("\ndone copy %d",hpp.reconMode);
	//Initialise output memory
	int sizeIm;
	mwSize dim[5] = { hpp.imageHeight,hpp.imageWidth,hpp.numChannel,hpp.numAngle,hpp.numFrame};
		
	switch (hpp.reconMode) {
	case 0:
		//HRI-Images
		DAS_beamforming_GPU_HRI(dpp, hpp);

		sizeIm = hpp.imageHeight*hpp.imageWidth*hpp.numFrame;
		dim[2] = hpp.numFrame; dim[3] = 1; dim[4] = 1;
		
		break;
	case 1:
		//LRI-Images
//         mexPrintf("Beamforming start");
		DAS_beamforming_GPU_LRI(dpp,hpp);
		
		sizeIm = hpp.imageHeight*hpp.imageWidth*hpp.numAngle*hpp.numFrame;
		dim[2] =hpp.numAngle ; dim[3] = hpp.numFrame; dim[4] = 1;
		
		break;
	case 2:
		//Channel images
		DAS_beamforming_GPU_CRI(dpp, hpp);
		sizeIm = hpp.imageHeight*hpp.imageWidth*hpp.numChannel*hpp.numAngle*hpp.numFrame;

		break;
	}

	//Return output images
//     mexPrintf("Transfering data");
	plhs[OUTPUTIMAGE] = mxCreateNumericArray(5, dim, mxSINGLE_CLASS, mxCOMPLEX); //create nd single 
	hpp.LRI_r = (float *)mxGetData(plhs[OUTPUTIMAGE]);
	hpp.LRI_i = (float *)mxGetImagData(plhs[OUTPUTIMAGE]);

//     mexPrintf("\n%d, %d,  %d, %d, %d", hpp.imageWidth, hpp.imageHeight, hpp.numChannel,hpp.numAngle,hpp.numFrame);
	for (i = 0; i<sizeIm; i++)
	{
		hpp.LRI_r[i] = dpp.h_LRI[i].x;
		hpp.LRI_i[i] = dpp.h_LRI[i].y;
	}
}
