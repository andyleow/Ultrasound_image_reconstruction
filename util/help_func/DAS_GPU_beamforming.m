% GPU beamforming wrapper developed by Chee Hau Leow
%  Property of ULIS group, Imperial College London
% Distribution of the program outside ULIS group are strictly prohibited
% Requests for permission should be addressed to c.leow12@ic.ac.uk/mengxing.tang@ic.ac.uk 

function[ImgData]= DAS_GPU_beamforming(location, Trans,ImagParam, FiltParam, pixelMap, progress)  
%%
% This function read the rf data from the location and beamform the data
% (plane wave/diverging wave) according to the predefine parameters.
% 1) Input: 
%   a) location: string contain the file location of the raw RF data
%   b) Trans: transducer properties enter via genTranducer.m function
%   c) ImagParam: struct data contain the following field
%         i) probe: Tranducer name
%         ii) startFrame: first frame of the imaging sequence
%        iii) endFrame: last frame of the imaging sequence
%       iv) numFiringPerFrame:  Number of transmission used to form one high resolution Image
%        v) skipFrame: Number of high resolution images to skip
%        vi) degX_tx: A vector contain the transmit steering angle in x direction
%        vii) degY_tx: A vector contain the transmit steering angle in y direction
%        viii) sumPulse: Pulse summing 0:Bmode, 1:PI all-no summing 2: PI sum-sum2pulse
%        ix) returnType: 'HRI' or 'LRI' reconstruction
%        x) fs: sampling frequency
%        xi) c: speed of sound
%        xii) gain: digital gain
%        xiii) delay: delay that account for transmit delay, receive delay and lens correction [s]
%       xiv) numSample: number of receive sample point per firing
%       xv) numChannel: number of receive channel
%       xvi) numAcq: number of acquiring frames per DMA transfer frame
%       xvii) numFrame: number of DMA transfer frame
%       xviii) senseCutoff: element sensitivity cutoff (0-1);
%   d) FiltParam: Struct data 
%       i) type:  'all','FIRHigh','FIRLow','FIRBand'
%       ii) Fpass1: 1st cutoff frequency
%       iii) Fpass2:  2nd cutoff frequency
%        iv) tempMedFilt: Temporal medium filter, 0=off; 1=on
%   e)pixelMap: Struct data 
%       i) pixelMapX:  lateral position
%       ii) pixelMapY:  elevation position
%       ii) pixelMapZ: axial position
%   f)tempMedFilt: 0 or 1 to enable temporal median filter 
% 2) Output: LRI, HRI or CRI images

%Specification for the GPU function:
% 
% Preallocate GPU resources and initialise imaging parameters
% if (rf_data is in integer format)
%  Eg : reconImg =  cuDAS_int(initMode,rf_data,degX_tx,degY_tx,actApertureX,...
%               actApertureY,actApertureZ,FilterCoef, pixelMapX,pixelMapY,...
%               pixelMapZ,delay,fs, ftx,c,focus,senseCutoff,ReconMode,gpuID);
% (rf_data is in single precision format) 
% reconImg =  cuDAS_single(initMode,rf_data,degX_tx,degY_tx,actApertureX,...
%               actApertureY,actApertureZ,FilterCoef, pixelMapX,pixelMapY,...
%               pixelMapZ,delay,fs, ftx,c,focus,senseCutoff,ReconMode,gpuID);

% Execute reconstruction without the need for parameters
% Eg :  reconImg=cuDAS_int(initMode,rf_data);
%
% Terminate GPU recon and free all GPU resource
% Eg :  cuDAS_int(initMode);

% Input: 
%   1)initMode: 0:Execute and/or initialie memory; 1:(Re-)initialise memory, 2:clear memory
%   2)rf_data: Input rf, (sample, channel,frame)/ (sample,channel,angle,frame)
%   3) degX_tx : A vector containing the transmit steering angle in x direction (in radian)
%   4) degY_tx : A vector containing the transmit steering angle in y direction (in radian)
%   5) actApertureX:vector containing the aperture position in x direction(m)
%   6) actApertureY:vector containing the aperture position in y direction(m)
%   7) actApertureZ:A vector containing the aperture position in z direction(m)
%   8) FilterCoef: A vector containing FIR filter coefficient 
%   9) pixelMapX: A vector containing the lateral image position (m)
%   10)pixelMapY: A vector containing the elevation image position (m)
%   11) pixelMapZ:A vector containing the axial image position (m)
%   12) delay: delay that account for transmit delay, receive delay and lens correction [s]
%   13) fs: sampling frequency
%   14) ftx: transmit frequency (for apodization function calculation)
%   15) c: speed of sound
%   16) focus: transmit focal point: 0 for plane wave, -ve for diverging wave 
%   17) senseCutoff: sensitivity cutoff (0-1)
%   17) ReconMode: 0 - LRI, 1 - HRI, 2 - Channel 
%   18) gpuID: gpu device number (if multiple exist)
% Output:
% reconImg: Reconstructed image (only if initMode=0) 

%%
% Define filtering coeff based on pre-filtering option
FiltParam.coef = computeFIR(ImagParam.fs,FiltParam);

%% load Receive Data
load(location,'RcvData','RcvLastFrame');
if isinteger(RcvData);
    RcvData=deserialize(RcvData);
end

if iscell(RcvData)
    sizeRcvData = size(RcvData{1});
else
    sizeRcvData=size(RcvData);
end

if length(sizeRcvData)<3
    sizeRcvData(3)=1;
end

try
    load(location,'splitPart');
end

if exist('splitPart','var')
    if splitPart~=1
        for i=2:splitPart
            eval(['load(location,''RcvData',num2str(i),'''',')']);
            eval(['RcvData',num2str(i),'=deserialize(RcvData',num2str(i),');']);
            eval([' RcvData{1}=cat(3,RcvData{1},','RcvData',num2str(i),'{1});']);
             eval(['clear RcvData',num2str(i)]);
        end
    end
elseif (sizeRcvData(3)<ImagParam.numFrame)
    if sizeRcvData(3)~=ImagParam.numFiringPerFrame
        load(location,'RcvData2');
        if isinteger(RcvData2);
            RcvData2=deserialize(RcvData2);
        end
        RcvData{1}=cat(3,RcvData{1},RcvData2{1});
        clear RcvData2;
    end
end

if iscell(RcvData)
    RcvData=RcvData{1};
end

if ~exist('RcvLastFrame','var')
    RcvLastFrame=sizeRcvData(3);
end
if RcvLastFrame<sizeRcvData(3)
   RcvData=circshift(RcvData,[0,0,sizeRcvData(3)-RcvLastFrame]);
end

%%  Restructure receive data

totalFrame = ceil((ImagParam.endFrame-ImagParam.startFrame+1)/ImagParam.skipFrame);
if totalFrame>ImagParam.numAcq*ImagParam.numFrame
    totalFrame=ImagParam.numAcq*ImagParam.numFrame;
end

if ImagParam.sumPulse>0
    n=ImagParam.numSample;
    k1=ImagParam.numFiringPerFrame*ImagParam.numAcq;
    switch ImagParam.sumPulse
        case 1
            disp('PI Pulses');
            ImagParam.numAcq=2*ImagParam.numAcq;
             totalFrame=ImagParam.numAcq*ImagParam.numFrame;
        case 2
            disp('Summing 2 Pulses');
            RcvData=RcvData(1:n*k1,:,:)+RcvData(n*k1+1:2*n*k1,:,:);
        case 3
            disp('Summing 3 Pulses ');
            RcvData=RcvData(1:n*k1,:,:)+RcvData(n*k1+1:2*n*k1,:,:)+RcvData(2*n*k1+1:3*n*k1,:,:);
        otherwise
            disp('Error');
            return;
    end
end

%% Temporal median filtering 
if FiltParam.tempMedFilt
    RcvData=permute(single(RcvData),[2,1,3]);
    nc=size(RcvData,1);
    nr=ImagParam.numSample;
    RcvData=reshape(RcvData,nc*nr*UserSet.na,[]);
    
    for i=2:size(RcvData,2)-1
%             disp(num2str(i))
            sliceInd= i-1:i+1;
            RcvData(:,i)=median(RcvData(:,sliceInd),2);
    end
       
    RcvData=reshape(RcvData,nc,[],UserSet.numFrame);
    RcvData=permute(RcvData,[2,1,3]);
    RcvData={int16(RcvData)};
    
%     temp=reshape(temp,nc,[],UserSet.numFrame);
%     temp=permute(temp,[2,1,3]);
%     clear RcvData;
%     RcvData=int16(temp);
%     clear temp
end

%% Define RF data location
Receive = repmat(struct('framenum', 1, ...
    'acqNum', 1, ...
    'startSample',1, ...
    'endSample',ImagParam.numSample),ImagParam.numFiringPerFrame*ImagParam.numAcq*ImagParam.numFrame,1);

if (ImagParam.numAcq==1 && ImagParam.numSample==sizeRcvData(1))
        for i = 1:ImagParam.numFrame*ImagParam.numFiringPerFrame % 3 acquisitions per frame
            Receive(i).framenum = i;
        end
else
    for i = 1:ImagParam.numFrame % 3 acquisitions per frame
        k = ImagParam.numFiringPerFrame*(i-1)*ImagParam.numAcq;
        for j = 1:ImagParam.numFiringPerFrame*ImagParam.numAcq
            Receive(k+j).framenum = i;
            Receive(k+j).acqNum = j;
            Receive(k+j).startSample=(j-1)*ImagParam.numSample+1;
            Receive(k+j).endSample=j*ImagParam.numSample;
        end
    end
end
%% Define Aperture Index and position
numChannel=sizeRcvData(2);
Channels=1:numChannel;

if numChannel<Trans.numElement
    actApertureX=Trans.position(Channels+Trans.startAperture-1,1);
    Trans.connector = Trans.connector(Channels+Trans.startAperture-1);
%     Channels=circshift(Channels,[0,-Trans.startAperture+1]);
else
     actApertureX=Trans.position(:,1);
     numChannel=Trans.numElement;    
end

actApertureY=Trans.position(:,2);
actApertureZ=Trans.position(:,3);

%%
% Generate pixel map of the imaging view
pixelMapX = pixelMap.upperLeft(1):pixelMap.dx:pixelMap.bottomRight(1);
pixelMapY = pixelMap.upperLeft(2):pixelMap.dy:pixelMap.bottomRight(2);
if isempty(pixelMapY)
    pixelMapY=0;
end
pixelMapZ = (pixelMap.upperLeft(3):pixelMap.dz:pixelMap.bottomRight(3))';
imageSize = [length(pixelMapZ),length(pixelMapX)] ;			% 5mm = 33pixel pixel.z

 if ~isfield(pixelMap,'focus')
     pixelMap.focus=0;
 end
 
 %%
% Initialise output image
switch (ImagParam.returnType)
    case 'LRI'
        ReconMode=1;
        ImgData = zeros(imageSize(1),imageSize(2),ImagParam.numFiringPerFrame,floor(totalFrame/ImagParam.skipFrame),'single');
    case 'CRI'
        ReconMode=2;
        ImgData = zeros(imageSize(1),imageSize(2),ImagParam.numChannel,ImagParam.numFiringPerFrame,floor(totalFrame/ImagParam.skipFrame),'single');
    otherwise
        ReconMode=0;
        ImgData = zeros(imageSize(1),imageSize(2),floor(totalFrame/ImagParam.skipFrame),'single');
end

% Turn on progress bar
if (strcmp(progress,'on'))
	h = waitbar(0,'Processing Frame #');
end

if isinteger(RcvData)
    rf_data=zeros(ImagParam.numSample,numChannel,ImagParam.numFiringPerFrame,'int16');
else
    rf_data=zeros(ImagParam.numSample,numChannel,ImagParam.numFiringPerFrame,'single');
end
for i = 1:totalFrame
    RcvIndex=(i-1)*ImagParam.skipFrame+ImagParam.startFrame;
    if (strcmp(progress,'on'))
        waitbar((i-1)/totalFrame,h,sprintf('Processing Plane Wave Frame #%d',i));
    else
        disp(['Processing Plane Wave Frame #',num2str(i)]);
    end
    
    %   % Update here to fit your own file name
    for j=1:ImagParam.numFiringPerFrame
        rf_data(:,:,j)=RcvData(Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).startSample:Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).endSample,...
            Trans.connector,Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).framenum);
    end
    
    if (i==1)
        if isinteger(rf_data)
            cuDAS_int(1,rf_data,ImagParam.degX_tx,ImagParam.degY_tx, actApertureX,actApertureY,actApertureZ,FiltParam.coef, pixelMapX,...
                pixelMapY,pixelMapZ,ImagParam.delay,ImagParam.fs, Trans.frequency,ImagParam.c, pixelMap.focus,ImagParam.senseCutoff,ReconMode,0);
        else
           cuDAS_single(1,rf_data,ImagParam.degX_tx,ImagParam.degY_tx, actApertureX,actApertureY,actApertureZ,FiltParam.coef, pixelMapX,...
                pixelMapY,pixelMapZ,ImagParam.delay,ImagParam.fs, Trans.frequency,ImagParam.c, pixelMap.focus,ImagParam.senseCutoff,ReconMode,0);
        end
    end
    
    if isinteger(rf_data)
        image1=cuDAS_int(0,rf_data);
    else
        image1=cuDAS_single(0,rf_data);
    end
    
    
    switch (ImagParam.returnType)
        case 'LRI'
            ImgData(:,:,:,i) = image1;
        case 'CRI'
            ImgData(:,:,:,:,i) = image1;
        otherwise
            ImgData(:,:,i) = image1;
    end
end
if isinteger(rf_data)
    cuDAS_int(-1);
else
    cuDAS_single(-1);
end
% Close progress bar if turned on
if (strcmp(progress,'on'))
    close (h);
end
end
