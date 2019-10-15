% GPU beamforming wrapper developed by Chee Hau Leow
%  Property of ULIS group, Imperial College London
% Distribution of the program outside ULIS group are strictly prohibited
% Requests for permission should be addressed to c.leow12@ic.ac.uk/mengxing.tang@ic.ac.uk 

function[ImgData]= Fourier_GPU_beamforming(location, Trans,ImagParam, FiltParam, pixelMap, progress)  
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
%        vi) deg_tx: A vector containing the transmit steering angle
%        vii) sumPulse: Pulse summing 0:Bmode, 1:PI all-no summing 2: PI sum-sum2pulse
%        viii) returnType: 'HRI' or 'LRI' reconstruction
%        ix) fs: sampling frequency
%         x) c: speed of sound
%         xi) gain: digital gain
%        xii) delay: delay that account for transmit delay, receive delay and lens correction [s]
%       xiii) numSample: number of receive sample point per firing
%       xiv) numChannel: number of receive channel
%       xv) numAcq: number of acquiring frames per DMA transfer frame
%       xvi) numFrame: number of DMA transfer frame
%   d) FiltParam: Struct data 
%       i) type:  'all','FIRHigh','FIRLow','FIRBand'
%       ii) Fpass1: 1st cutoff frequency
%       iii) Fpass2:  2nd cutoff frequency
%        iv) tempMedFilt: Temporal medium filter, 0=off; 1=on
%   e)pixelMap: Struct data 
%       i) pixelMapX:  lateral position
%       ii) pixelMapZ: axial position
%   f)tempMedFilt: 0 or 1 to enable temporal median filter 
% 2) Output: LRI, HRI or CRI images

%Specification for the GPU function: 
%  Eg : reconImg = cuDAS_mainV1(initMode,rf_data, deg_tx, actAperture,FiltCoef, 
% pixelMapX, t0, pixelMapZ,fs,ftx,c,focus,ReconMode,gpuID);

% Input: 
%   1)initMode: 0:Execute and/or initialie memory; 1:(Re-)initialise memory, 2:clear memory
%   2)rf_data: Input rf, multiPW(sample, channel,angle,frame)/singlePW(sample, channel, frame)
%   3) deg_tx : A vector containing the transmit steering angle (in radian)
%   4) actAperture:A vector containing the aperture position
%   5) FilterCoef: A vector containing FIR filter coefficient
%   6) pixelMapX: A vector containing the lateral image position
%   7) pixelMapZ:A vector containing the axial image position
%   8) delay: delay that account for transmit delay, receive delay and lens correction [s]
%   9) fs: sampling frequency
%   10) ftx: transmit frequency (for apodization function calculation)
%   11) c: speed of sound
%   12) focus: transmit focal point: 0 for plane wave, -ve for diverging wave 
%   13) ReconMode: 0 for LRI and 1 for HRI
%   14) gpuID: gpu device number (if multiple exist)
% Output:
% reconImg: Reconstructed image
%

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
    k1=ImagParam.NumFiringPerFrame*ImagParam.numAcq;
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
    RcvData=permute(single(RcvData{1}),[2,1,3]);
    nc=size(RcvData,1);
    nr=Intermediate.numRcvSamples;
    RcvData=reshape(RcvData,nc*nr*UserSet.na,[]);
    
    for i=2:size(RcvData,2)-1
%             disp(num2str(i))
            sliceInd= i-1:i+1;
            RcvData(:,i)=median(RcvData(:,sliceInd),2);
    end
       
    RcvData=reshape(RcvData,nc,[],UserSet.NumFrame);
    RcvData=permute(RcvData,[2,1,3]);
    RcvData={int16(RcvData)};
    
%     temp=reshape(temp,nc,[],UserSet.NumFrame);
%     temp=permute(temp,[2,1,3]);
%     clear RcvData;
%     RcvData{1}=int16(temp);
%     clear temp
end

%% Define RF data location

Receive = repmat(struct('framenum', 1, ...
    'acqNum', 1, ...
    'startSample',1, ...
    'endSample',ImagParam.numSample),ImagParam.numFiringPerFrame*ImagParam.numAcq*ImagParam.numFrame,1);

if (ImagParam.numAcq==1)
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
ChannelFinal=Channels;
if numChannel<Trans.numElement
    actAperture=Trans.position(Channels+Trans.startAperture-1);
    Channels=circshift(Channels,[0,-Trans.startAperture+1]);
else
     actAperture=Trans.position(:,1)';
     if isfield(Trans,'connector')
        [Channels,ChannelFinal] = sort(Trans.connector);
        numChannel=Trans.numElement;    
     end
end

%%
% Generate pixel map of the imaging view
pixelMapX = pixelMap.upperLeft(1):pixelMap.dx:pixelMap.bottomRight(1);
pixelMapZ = (pixelMap.upperLeft(2):pixelMap.dz:pixelMap.bottomRight(2))';
imageSize = [length(pixelMapZ),length(pixelMapX)] ;			% 5mm = 33pixel pixel.z

 if ~isfield(pixelMap,'focus')
     pixelMap.focus=0;
 end
 
 
 %%
% Initialise output image
% imageHeight=2^nextpow2(pixelMap.height);
% imageWidth = 2^(nextpow2(pixelMap.width)-1);
switch (ImagParam.returnType)
    case 'LRI'
        ReconMode=1;
        ImgData = zeros(imageSize(1),imageSize(2),ImagParam.numFiringPerFrame,floor(totalFrame/ImagParam.skipFrame),'single');
    otherwise
        ReconMode=0;
         ImgData = zeros(imageSize(1),imageSize(2),floor(totalFrame/ImagParam.skipFrame),'single');
end

if pixelMap.focus~=0
    error('Fourier beamforming not supporting diverging/focus reconstruction')
    return 
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
        rf_data(:,ChannelFinal,j)=RcvData(Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).startSample:Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).endSample,...
            Channels,Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).framenum);
    end
    
      if (i==1)
        if isinteger(rf_data)
            cuFBF_int(1,rf_data,ImagParam.deg_tx, actAperture,FiltParam.coef, pixelMapX,...
                pixelMapZ,ImagParam.delay,ImagParam.fs, Trans.frequency,ImagParam.c, pixelMap.focus,ReconMode,0);
        else
            cuFBF_single(1,rf_data,ImagParam.deg_tx, actAperture,FiltParam.coef, pixelMapX,...
                pixelMapZ,ImagParam.delay,double(ImagParam.fs), Trans.frequency,ImagParam.c, pixelMap.focus,ReconMode,0);
        end
    end
        
   if isinteger(rf_data)
        image1=cuFBF_int(0,rf_data);
    else
        image1=cuFBF_single(0,rf_data);
    end
    
    switch (ImagParam.returnType)
        case 'LRI'
            ImgData(:,:,:,i) = image1;
        otherwise
            ImgData(:,:,i) = image1;
    end
end

if isinteger(rf_data)
    cuFBF_int(-1);
else
    cuFBF_single(-1);
end

% Close progress bar if turned on
if (strcmp(progress,'on'))
    close (h);
end
end
