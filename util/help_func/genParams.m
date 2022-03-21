function [UserSet,Trans,ImagParam] = genParams(varargin)

% Generate transducer and Imaging Parameters for reconstruction
% Input:
% a) UserSet: structure contain acquisition paramters
%       i) probe (Optional): transducer name,  if not exist, must specified
%           in probe parameter
%       ii)  numFrame/NumFrame:; number of DMA transfer frame
%       iii) numAcq: number of acquiring frames per DMA transfer frame
%       iv) na: number of angle in x direction
%       v) na2: number of angle in y direction (only exist in 3D imaging)
%       vi) angleRange: transmit angle Range in x direction
%       vii) angleRange2: transmit angle Range in y direction (only exist in 3D imaging)
%       viii) numChannel (Optional): number of receive channel
%       ix) numSample (Optional): number of receive sample per firing
%       x) sampleDepthMM: receive sample in (mm). Used to calculate
%           numSample if numSample is not specified
%       xi) txFreq/TXFreq: transmit frequency
%  b) Probe : transducer name,  if not exist in UserSet
%  c) Recon (Optional): Struct contain recon information
%       i) startFrame: first frame of the imaging sequence
%       ii) endFrame: last frame of the imaging sequence
%       iii) numFiringPerFrame:  Number of transmission used to form one high resolution Image
%       iv) skipFrame: Number of high resolution images to skip
%       vi) sumPulse: Pulse summing 0:Bmode, 1:PI all-no summing 2: PI sum-sum2pulse
%        viii) type: HRI' or 'LRI' reconstruction
%
%   Output:
%  d) UserSet: structure contain fully populated acquisition paramters
%  e) Trans: structure contain transducer parameters
%  f) ImagParam: struct data contain the following field
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

%%    
if nargin==1
    UserSet=varargin{1};
    Probe =[];
    Recon = [];
elseif nargin==2
    UserSet=varargin{1};
    Probe=varargin{2};
    Recon=[];
else
    UserSet=varargin{1};
    Probe=varargin{2};
     Recon=varargin{3};
end

%Default value
ImagParam = struct('probe',[],'startFrame',1,'endFrame',1,'numFiringPerFrame',1,'skipFrame',1,...
    'degX_tx',0,'degY_tx',0, 'sumPulse', 0,'returnType','HRI', 'fs',[],'c',1540,'gain',1,'delay',[], ...
    'numSample',[],'numChannel',128,'numAcq',1, 'numFrame',1);

% Transducer parameter
if isfield(UserSet, 'probe')
    Trans = genTransducer(UserSet.probe);
else
    Trans = genTransducer(Probe.name);
end
ImagParam.probe = Trans.name;

%standardise UserSet format
if isfield(UserSet,'NumFrame')
    UserSet.numFrame =UserSet.NumFrame;
    UserSet = rmfield(UserSet,'NumFrame');
end

if isfield(UserSet,'AngleRange')
    UserSet.angleRange =UserSet.AngleRange;
    UserSet = rmfield(UserSet,'AngleRange');
end

if isfield(UserSet,'TXFreq')
    UserSet.txFreq =UserSet.TXFreq;
    UserSet = rmfield(UserSet,'TXFreq');
end

% Imaging parameter
if isfield(Recon,'startFrame')
    ImagParam.startFrame =Recon.startFrame;
else
    ImagParam.startFrame = 1;
end

if isfield(Recon,'endFrame')
    ImagParam.endFrame = Recon.endFrame;
else
    ImagParam.endFrame = UserSet.numAcq*UserSet.numFrame;
end

if isfield(Recon,'numFiringPerFrame')
    ImagParam.numFiringPerFrame = Recon.numFiringPerFrame;
else
    ImagParam.numFiringPerFrame = UserSet.na;
end

if isfield(Recon,'skipFrame')
    ImagParam.skipFrame = Recon.skipFrame;  % 1 = no skip, 2 = skip 1 frame, etc
else
    ImagParam.skipFrame = 1; % 1 = no skip, 2 = skip 1 frame, etc
end

if isfield(Recon,'sumPulse')
    ImagParam.sumPulse=Recon.sumPulse;       %Pulse summing 0:Bmode, 1:PI all-no summing 2: PI sum-sum2pulse
else
    ImagParam.sumPulse=0;
end

if isfield(Recon,'type')
    ImagParam.returnType=Recon.type;       %Pulse summing 0:Bmode, 1:PI all-no summing 2: PI sum-sum2pulse
end
 
% Sampling frequency
if isfield(UserSet,'fs')
    ImagParam.fs = double(UserSet.fs);
else
    ImagParam.fs = Trans.frequency*UserSet.samplePerWave;
end

% Transmit angle
if UserSet.na>1
    if isfield(UserSet,'na2') && isfield(UserSet,'angleRange2')
        if UserSet.na2>1
            for i=1:UserSet.na
                for j=1:UserSet.na2
                    index =(i-1)*UserSet.na+j;
                    ImagParam.degX_tx(index)=(-UserSet.angleRange/2+(i-1)*(UserSet.angleRange)/(UserSet.na-1))*pi/180;
                    ImagParam.degY_tx(index)=(-UserSet.angleRange2/2+(j-1)*(UserSet.angleRange2)/(UserSet.na2-1))*pi/180;
                end
            end
        end
    else
        for i=1:UserSet.na
            ImagParam.degX_tx(i)=(-UserSet.angleRange/2+(i-1)*(UserSet.angleRange)/(UserSet.na-1))*pi/180;
            ImagParam.degY_tx(i)=0;
        end
    end
end

% %Receive angle (for doppler processing)
% ImagParam.deg_rx = ones(ImagParam.numFiringPerFrame,1)*0;

% Delay compensation
ImagParam.delay = (2*Trans.lensCorrection/ImagParam.c)*ones(ImagParam.numFiringPerFrame,1); %zeros(ImagParam.NumFiringPerFrame,1);

%Receive channel
if isfield(UserSet,'numChannel')
    ImagParam.numChannel = UserSet.numChannel;
else
    if Trans.name(1)=='M'
        ImagParam.numChannel = Trans.numElement;
    end
end

% Automatic find the first channel if aperture greater than receive channel
% or reduce ImagParam.numChannel to map with the transducer channel
Trans.startAperture =1; 
if Trans.numElement>ImagParam.numChannel
    Trans.startAperture=(Trans.numElement-ImagParam.numChannel)/2+1;
elseif Trans.numElement<ImagParam.numChannel
    ImagParam.numChannel = Trans.numElement;
end


%Calculate the number of sample per firing is not specified
if ~isfield(UserSet,'numSample')
    switch Trans.name(1)
        case 'C' 
            Intermediate.sampleDepthWL=ceil(UserSet.sampleDepthMM*1e-3*(Trans.frequency/1540));
            Intermediate.maxAcqLength = sqrt(Intermediate.sampleDepthWL^2 + (ImagParam.numChannel*Trans.pitch*(Trans.frequency/ImagParam.c))^2);
            Intermediate.wlsPer128 = ImagParam.numChannel/(UserSet.samplePerWave*2); % wavelengths in 128 samples for PI
            Intermediate.numRcvSamples = 2*(Intermediate.wlsPer128*ceil(Intermediate.maxAcqLength/Intermediate.wlsPer128))*UserSet.samplePerWave;
        case 'L' 
            Intermediate.sampleDepthWL=ceil(UserSet.sampleDepthMM*1e-3*(Trans.frequency/1540));
            Intermediate.maxAcqLength = sqrt(Intermediate.sampleDepthWL^2 + (ImagParam.numChannel*Trans.pitch*(Trans.frequency/ImagParam.c))^2);
            Intermediate.wlsPer128 = ImagParam.numChannel/(UserSet.samplePerWave*2); % wavelengths in 128 samples for PI
            Intermediate.numRcvSamples = 2*(Intermediate.wlsPer128*ceil(Intermediate.maxAcqLength/Intermediate.wlsPer128))*UserSet.samplePerWave;
        case 'P'
            Intermediate.sampleDepthWL=ceil(UserSet.sampleDepthMM*1e-3*(Trans.frequency/1540));
            Intermediate.maxAcqLength = sqrt(Intermediate.sampleDepthWL^2 + (Trans.numElement*Trans.pitch*(Trans.frequency/ImagParam.c))^2+2* Intermediate.sampleDepthWL*cos( UserSet.scanRangle-pi/2));
            Intermediate.wlsPer128 = ImagParam.numChannel/(UserSet.samplePerWave*2); % wavelengths in 128 samples for PI
            Intermediate.numRcvSamples = 2*(Intermediate.wlsPer128*ceil(Intermediate.maxAcqLength/Intermediate.wlsPer128))*UserSet.samplePerWave;
        case 'M'  
            Intermediate.sampleDepthWL=ceil(UserSet.sampleDepthMM*1e-3*(Trans.frequency/1540));
            Intermediate.maxAcqLength = sqrt(Intermediate.sampleDepthWL^2 + (2*max(Trans.position(:))*(Trans.frequency/ImagParam.c))^2);
            Intermediate.wlsPer128 = 256/(UserSet.samplePerWave*2); % wavelengths in 128 samples for PI
            Intermediate.numRcvSamples = 2*(Intermediate.wlsPer128*ceil(Intermediate.maxAcqLength/Intermediate.wlsPer128))*UserSet.samplePerWave;
    end
    UserSet.numSample= Intermediate.numRcvSamples;
    ImagParam.numSample=Intermediate.numRcvSamples;
else
    ImagParam.numSample= UserSet.numSample;
end

ImagParam.numAcq = UserSet.numAcq;
ImagParam.numFrame = UserSet.numFrame;

return



