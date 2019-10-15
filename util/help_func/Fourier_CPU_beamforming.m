% GPU beamforming wrapper developed by Chee Hau Leow
%  Property of ULIS group, Imperial College London
% Distribution of the program outside ULIS group are strictly prohibited
% Requests for permission should be addressed to c.leow12@ic.ac.uk/mengxing.tang@ic.ac.uk 

function[ImgData]= Fourier_CPU_beamforming(location, Trans,ImagParam, FiltParam, pixelMap, progress)  
%%
% This function read the rf data from the location and beamform the data
% (plane wave only) according to the predefine parameters.
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
%       xvii) senseCutoff: element sensitivity cutoff (0-1);
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
%  Eg : [BF_img] = FBF_CPU_Func(rf_data, tx_angle, actElePosX,filterCoef, pixelMapX, 
%       pixelMapZ,t0,fs,ftx,c,focus,ReconMode)

% Input:
%   1)rf_data: Input rf, multiPW(sample, channel,angle,frame)/singlePW(sample, channel, frame)
%   2) tx_angle : A vector containing the transmit steering angle (in radian)
%   3) actElePosX:A vector containing the aperture position
%   4) filterCoef: A vector containing FIR filter coefficient
%   5) pixelMapX: A vector containing the lateral image position
%   6) pixelMapZ:A vector containing the axial image position
%   7) t0: delay that account for transmit delay, receive delay and lens correction [s]
%   8) fs: sampling frequency
%   9) ftx: transmit frequency (for apodization function calculation)
%   10)  c: speed of sound
%   11) focus: transmit focal point: 0 for plane wave, -ve for diverging wave
%   12)  ReconMode: 0 for LRI and 1 for HRI
% Output:
% BF_img: Reconstructed image

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

%%
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
% Generate pixel map of the imaging view
pixelMapX = pixelMap.upperLeft(1):pixelMap.dx:pixelMap.bottomRight(1);
pixelMapZ = (pixelMap.upperLeft(2):pixelMap.dz:pixelMap.bottomRight(2))';
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
        rf_data(:,ChannelFinal,j)=RcvData(Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).startSample:Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).endSample,...
            Channels,Receive((RcvIndex-1)*ImagParam.numFiringPerFrame+j).framenum);
    end
    
    image1=FBF_CPU_Func(rf_data,ImagParam.deg_tx, actAperture,FiltParam.coef, pixelMapX,...
            pixelMapZ,ImagParam.delay,ImagParam.fs, Trans.frequency,ImagParam.c, pixelMap.focus,ReconMode);
    
        
    switch (ImagParam.returnType)
        case 'LRI'
            ImgData(:,:,:,i) = image1;
        otherwise
            ImgData(:,:,i) = image1;
    end
end

% Close progress bar if turned on
if (strcmp(progress,'on'))
    close (h);
end
end
