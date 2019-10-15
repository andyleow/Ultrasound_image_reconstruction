%Script written by Chee Hau Leow
%email: c.leow12@ic.ac.uk
%GPU beamforming adapted in external function

%% Include path
%cuda, Serialisation and savefast function path
cd('/home/chleow/Matlab/Vantage-4.0.1-1903121200');
activate; 

%Serialisation and savefast function path
% addpath(genpath('C:\Users\chl912\Documents\Matlab_ExFun\02 Matlab Savefast\savefast'));
% addpath(genpath('C:\Users\chl912\Documents\VerasonicsSoftware\CustomScript\prototype\L22-14v\gpu_prototype\cuda'));
addpath(genpath('/home/chleow/Matlab/externFunc'));
addpath(genpath('/home/chleow/Matlab/Vantage-4.0.1-1903121200/myFunc/cuda'));

%%
close all;
clear;
% clearvars -except RcvData; %Uncomment if sim=2;

savedir=['/home/chleow/Matlab/Vantage-4.0.1-1903121200/TempData\',datestr(now,'yymmdd'),'\'];
mkdir(savedir);

if(exist('MatFiles/Preset.mat','file'))
    load('MatFiles/Preset.mat');
else
   Preset.Voltage=10;
    Preset.Range=60;
    Preset.SenseCut=0.5;
    Preset.CompressMethod='log'; %or 'power'
    Preset.Compress=50;
    Preset.DG=2;
    Preset.TGC=[613,660,760,850,945,1023,1023,1023];   
end

UserSet.sampleDepthMM=20;   %[15,20,25,30,35,40,45,50] in mm
UserSet.TOF=40;             %[54,58,62,66,70,78,]%[51,51,57,63,70,77,83,89]in us     
UserSet.FR=600;        %FrameRate
UserSet.numAcq=9; 
UserSet.NumFrame=40;            %Optimum Total acq=1000-1800 acq/ need to be even number
UserSet.Sim=0;                  %0=hardware, 1=sim, 2=data-rerun
UserSet.SimDownSample=1;

UserSet.TXFreq=18;       %Transmit Frequency(MHz)
UserSet.TXFocus= 0;
UserSet.numCycle=1;     %NumHalfCycle
UserSet.na=11;           %Num Compounded PW
UserSet.AngleRange=10;  %Max Range 
UserSet.samplePerWave=4;
UserSet.SVDEn=0;
UserSet.ReconSkipFrame=1;

%%
% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;            % warning level 0:3(error,warning,status, debug)
Resource.Parameters.initializeOnly = 0;     %set to initialise parameter w/0 running the hardware
Resource.Parameters.simulateMode = UserSet.Sim;

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L12-5 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 20; 

%% set dtheta to range over +/- 18 degrees.
if (UserSet.na > 1)
    dtheta = (UserSet.AngleRange*pi/180)/(UserSet.na-1); 
    startAngle=-UserSet.AngleRange*pi/180/2;
else
    dtheta = 0;
    startAngle=0; 
end

%%  Calculate acq depth,TOF, PRF
sampleDepthWL=ceil(UserSet.sampleDepthMM*(Trans.frequency/1.54));
Intermediate.maxAcqLength = sqrt(sampleDepthWL^2 + (Resource.Parameters.numRcvChannels*Trans.spacing)^2);
Intermediate.wlsPer128 = Resource.Parameters.numRcvChannels/(UserSet.samplePerWave*2); % wavelengths in 128 samples for PI
Intermediate.numRcvSamples = 2*(Intermediate.wlsPer128*ceil(Intermediate.maxAcqLength/Intermediate.wlsPer128))*UserSet.samplePerWave;
Intermediate.maxAcqLength=Intermediate.numRcvSamples/(UserSet.samplePerWave*2);  % final acquisition depth in wavelength

% Check PRF
Intermediate.TOF=ceil(2*Intermediate.maxAcqLength/Trans.frequency*1.1);

if UserSet.TOF<Intermediate.TOF
    UserSet.TOF=Intermediate.TOF;
end

if (1/UserSet.FR)*1e6 > (UserSet.TOF*UserSet.na) 
    PRF=round((1/UserSet.FR)*1e6-(UserSet.TOF*UserSet.na));
else
    PRF=0;
    UserSet.FR= 1/(UserSet.TOF*UserSet.na)*1e6;
end

%%
% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.125]; %x,y,z
PData.Size(1) = 2^nextpow2((sampleDepthWL))/PData.PDelta(3); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2)=((Resource.Parameters.numRcvChannels*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Resource.Parameters.numRcvChannels-1)/2,0,0]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation=-0.5;
Media.function = 'movePoints';

%%
% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = Intermediate.numRcvSamples*UserSet.numAcq*UserSet.na; %Specify for max range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = UserSet.NumFrame;        % 40 frames used for RF cineloop.

% Resource.InterBuffer(1).datatype = 'complex single';
% Resource.InterBuffer(1).pagesPerFrame=UserSet.numAcq;
% Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.

Resource.ImageBuffer(1).datatype = 'double';
if UserSet.Sim<2,Resource.ImageBuffer(1).numFrames = 10; 
else Resource.ImageBuffer(1).numFrames=fix(UserSet.numAcq*UserSet.NumFrame/UserSet.SimDownSample); end

Resource.DisplayWindow(1).Title = ['HFR_FlashAngles FrameRate=', num2str(UserSet.FR)];
Resource.DisplayWindow(1).pdelta = 1;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 10;   %Cineloop frame
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

%%
% Specify Transmit waveform structure. 
TW(1).type = 'parametric';
TW(1).Parameters = [UserSet.TXFreq,0.67,UserSet.numCycle,1];   % A, B, C, D  for 8.1 MHz transmit frequency

% Specify TX structure array.  
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ... 
                   'Apod', ones(1,Resource.Parameters.numTransmit), ...
                   'focus', UserSet.TXFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, UserSet.na);
               
% - Set event specific TX attributes.
for n = 1:UserSet.na   % 3*na transmit events
    TX(n).Steer = [(startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n)); % use 'TransmitOnAllElements' flag 
end

%%
% Specify TGC Waveform structure.
TGC.CntrlPts = Preset.TGC;
TGC.rangeMax = sampleDepthWL;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
Receive = repmat(struct('Apod', ones(1,Resource.Parameters.numRcvChannels), ...
                        'startDepth', 1, ...
                        'endDepth', Intermediate.maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,UserSet.na*UserSet.numAcq*UserSet.NumFrame);
                    
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames  % 3 acquisitions per frame
    k = UserSet.na*(i-1)*UserSet.numAcq;
%      Receive(k+1).callMediaFunc = 1;
    for j = 1:UserSet.na*UserSet.numAcq
        if mod(j,UserSet.na)==0
            Receive(k+j).callMediaFunc = 1;
        end
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

%%
% Specify Process structure array.
Process(1).classname = 'External';
Process(1).method = 'GPU_init';
Process(1).Parameters = {'srcbuffer','none'}; % name of buffer to process.
           

Process(2).classname = 'External';
Process(2).method = 'GPU_beamforming';
Process(2).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','image',...
                         'dstbufnum',1,...
                         'dstframenum',-2};
                     
pers = 0;
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',Preset.DG,...            % pgain is image processing gain
                         'reject',0,...
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','log',...
                         'compressFactor',Preset.Compress,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};                  

%%                     
% Specify SeqControl structure arrays.
SeqControl(1).command = 'timeToNextAcq';  % time between pulse
SeqControl(1).argument = UserSet.TOF;  % 20 usec
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = PRF+UserSet.TOF;  % 20 msec
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'jump'; % jump back to start.
SeqControl(4).argument = 2;
nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events
Event(n).info = 'Initialise GPU';
Event(n).tx = 0;         % use 1st TX structure.
Event(n).rcv = 0;    % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 1;    % no processing
Event(n).seqControl = 0; % time between syn. aper. acqs.
n = n+1;

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    l=(i-1)*UserSet.numAcq*UserSet.na;
    for h=1:UserSet.numAcq
        k = UserSet.na*(h-1);
        for j = 1:UserSet.na
            Event(n).info = 'tx-rx';
            Event(n).tx = j;         % use 1st TX structure.
            Event(n).rcv = l+k+j;    % use 1st Rcv structure.
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 1; % time between syn. aper. acqs.
            n = n+1;
        end
        Event(n-1).seqControl = 2; % Delay between frame
    end
    
    Event(n-1).seqControl = [2,nsc]; % use SeqControl structs defined below.
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;            %Kind of redundant(To be tested)
   
       if UserSet.Sim<2
           Event(n).info = 'GPU Reconstruction';
           Event(n).tx = 0;         % no transmit
           Event(n).rcv = 0;        % no rcv
           Event(n).recon = 0;      % reconstruction
           Event(n).process = 2;    % processing
           Event(n).seqControl=[3,nsc,nsc+1];
%            if floor(i/5) == i/5  
%        Event(n).seqControl = 3;
            SeqControl(nsc).command='waitForTransferComplete';
            SeqControl(nsc+1).command='markTransferProcessed';
            SeqControl(nsc).argument=nsc-1;
            SeqControl(nsc+1).argument=nsc-1;
            nsc = nsc + 2;
%            end
           n = n+1;
           
           Event(n).info = 'Process Image';
           Event(n).tx = 0;         % no transmit
           Event(n).rcv = 0;        % no rcv
           Event(n).recon = 0;      % reconstruction
           Event(n).process = 3;    % processing
           Event(n).seqControl = 3;
           n=n+1;         
       else
           for j= 1:fix(UserSet.numAcq/UserSet.SimDownSample)
               Event(n).info = 'Sim Reconstruct';
               Event(n).tx = 0;         % no transmit
               Event(n).rcv = 0;        % no rcv
               Event(n).recon = 0;      % reconstruction
               Event(n).process = 1;    % processing
               Event(n).seqControl = 3;
               n = n+1;
           end
       end
end

if UserSet.Sim<2
    Event(n).info = 'Jump back to first event';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0;
    Event(n).seqControl = 4; % jump command
end

%%
freeze=0;

UI(1).Control = {'UserB2','Style','VsPushButton','Label','ReconAll'};
UI(1).Callback = text2cell('%ReconRFCallBack');

% Preset Transmit Voltage
UI(2).Statement ='[result,hv] = setTpcProfileHighVoltage(evalin(''base'',''Preset.Voltage'')'',1);';
UI(3).Statement ='hv1Sldr = findobj(''Tag'',''hv1Sldr'');';
UI(4).Statement ='set(hv1Sldr,''Value'',hv);';
UI(5).Statement ='hv1Value = findobj(''Tag'',''hv1Value'');';
UI(6).Statement ='set(hv1Value,''String'',num2str(hv,''%.1f''));';

UI(7).Control = {'UserB1','Style','VsPushButton','Label','saveRF'};
UI(7).Callback = text2cell('%saveRFCallback');

%% How many External Functions?
NumOfFunc = 2;
for No = 1:NumOfFunc
    EF(No).Function = text2cell(['%EF#',num2str(No)]);
end

%%
% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = UserSet.numAcq/UserSet.na;

% Save all the structures to a .mat file.
save('MatFiles/L12-3v_FlashAngles_HFR');
filename='L12-3v_FlashAngles_HFR'; VSX;
return

%%
%ReconRFCallBack
if evalin('base','freeze')==0   % no action if not in freeze
    return
end

Control.Command = 'copyBuffers'; % copyBuffers does a sequenceStop
runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.
RcvData=RcvData{1};

GPUParams = evalin('base', 'GPUParams');
UserSet=evalin('base','UserSet');
GPUParams.numAcq = UserSet.numAcq;
GPUParams.numFrame=UserSet.NumFrame*UserSet.numAcq;
GPUParams.reconFrames=1:GPUParams.numFrame;
GPUParams.rf_data=zeros(GPUParams.numSample,GPUParams.numChannel,GPUParams.na,GPUParams.numFrame,'int16');

%Initialise GPU memory      
cuFBF_int(1,GPUParams.rf_data(:,:,:,1),GPUParams.deg_tx, GPUParams.actAperture,GPUParams.filterCoef, GPUParams.pixelMapX,...
            GPUParams.pixelMapZ,GPUParams.delay,GPUParams.fs, GPUParams.ftx,GPUParams.c, GPUParams.txFocus,...
            GPUParams.reconMode,GPUParams.gpuID);
RcvData=permute(reshape(RcvData,GPUParams.numSample,GPUParams.na,GPUParams.numAcq,[],UserSet.NumFrame),[1,4,2,3,5]);
GPUParams.rf_data=RcvData(:,GPUParams.channelInd,:,GPUParams.reconFrames);

h = waitbar(0,'Processing Frame #');
for i=1:GPUParams.numFrame
    if mod(i,10)==0
        waitbar((i-1)/GPUParams.numFrame,h,sprintf('Processing Frame %d / %d',i,GPUParams.numFrame));
    end
    GPUParams.IQData(:,:,i)=cuFBF_int(0,GPUParams.rf_data(:,:,:,i));
end
close (h);

% Display and save images
savedir=evalin('base','savedir');
vid1=VideoWriter([fullfile(savedir,datestr(now,'yymmdd_HHMMSS')),'_Bmode.avi']);
vid1.FrameRate = 100;
vid1.Quality = 75;

open(vid1);
f=figure

normMax=20*log10(max(abs(GPUParams.IQData(:))));

for i=1:GPUParams.numFrame
    imagesc(GPUParams.pixelMapX*1e3,GPUParams.pixelMapZ*1e3,20*log10(abs(GPUParams.IQData(:,:,i))),[normMax-50,normMax]);
    title(['Frame=',num2str(i)]);
    xlabel('mm');
    ylabel('mm');
    axis image;
    colormap(gray);
    
    F(i)=getframe(gcf);
    drawnow
    writeVideo(vid1,F(i));
end
%
close(vid1);
close(f);
return
%ReconRFCallBack

%saveRFCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
if simMode == 2
    return
end

if evalin('base','freeze')==0   % no action if not in freeze
    return
end

%Copy buffer back to RcvData;
Control.Command = 'copyBuffers'; % copyBuffers does a sequenceStop
runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.

%save RcvData
RcvData=evalin('base','RcvData(1)');

S=whos('RcvData');
splitPart=ceil(S.bytes/2^32);
% 
if splitPart~=1
    for i=2:splitPart
        part=ceil((i-1)*size(RcvData{1},3)/splitPart+1):ceil(i*size(RcvData{1},3)/splitPart);
        eval(['RcvData',num2str(i),'{1}','=RcvData{1}(:,:,part);']);
    end
    part=ceil((1-1)*size(RcvData{1},3)/splitPart+1):ceil(1*size(RcvData{1},3)/splitPart);
    RcvData{1}=RcvData{1}(:,:,part);
end

UserSet=evalin('base','UserSet');
%Find transmit voltage
hv1Sldr = findobj('Tag','hv1Sldr');
UserSet.voltage1=get(hv1Sldr,'Value');

fname = [datestr(now,'yymmdd_HHMMSS'),'PW_PI']
RcvLastFrame=size(RcvData{1},3);
if (~evalin('base','simButton'))
    RcvLastFrame=Resource.RcvBuffer(1).lastFrame;
end
direc=evalin('base','savedir');
RcvData=serialize(RcvData);

if splitPart~=1
    for i=2:splitPart
        eval(['RcvData',num2str(i),'=serialize(RcvData',num2str(i),');']);
    end
    
    string=['savefast([direc,fname],''RcvData'',''UserSet'',''RcvLastFrame'',''splitPart',''''];
    for i=2:splitPart
        var=['''RcvData',num2str(i),''''];
       string=[string,',',var];
    end
    string=[string,')'];
    eval(string);
else  
    savefast([direc,fname],'RcvData','UserSet','RcvLastFrame');
end

disp('.......done saving');
return
%saveRFCallback

%% External functions with process object

%EF#1
GPU_init()
% keyboard
GPUParams=struct('rf_data',[],'deg_tx',[], 'actAperture',[],'filterCoef',[],...
    'pixelMapX',[], 'pixelMapZ',[],'fs',[],'c',[],'delay',[],'txFocus',[],...
    'reconMode',0,'gpuID',0);

%Initialise GPUParams
Intermediate=evalin('base','Intermediate');
UserSet=evalin('base','UserSet');
Resource=evalin('base','Resource');
Trans=evalin('base','Trans');
PData= evalin('base','PData');

GPUParams.numSample=Intermediate.numRcvSamples;
GPUParams.numAcq=UserSet.numAcq;
GPUParams.na=UserSet.na;
GPUParams.numFrame=max(floor(UserSet.numAcq/UserSet.ReconSkipFrame),1);

if Trans.numelements<Resource.Parameters.numRcvChannels
    GPUParams.numChannel=Trans.numelements;
    GPUParams.channelInd = Trans.Connector;
else
    GPUParams.numChannel=Resource.Parameters.numRcvChannels;
    GPUParams.channelInd = Trans.Connector;
end

GPUParams.reconSkipFrame =UserSet.ReconSkipFrame;
GPUParams.reconFrames= 1:GPUParams.reconSkipFrame:GPUParams.numAcq;
GPUParams.rf_data=zeros(GPUParams.numSample,GPUParams.numChannel,GPUParams.na,GPUParams.numFrame,'int16');

GPUParams.deg_tx(1)=0;
if UserSet.na>1
    for i=1:UserSet.na
        GPUParams.deg_tx(i)=(-UserSet.AngleRange/2+(i-1)*(UserSet.AngleRange)/(UserSet.na-1))*pi/180;
    end
end

GPUParams.actAperture=Trans.ElementPos(:,1)*1e-3;
GPUParams.filterCoef= fir1(51,Trans.Bandwidth/(Trans.frequency*2));

GPUParams.fs=Trans.frequency*1e6*4;
GPUParams.ftx =UserSet.TXFreq*1e6;
GPUParams.c = Resource.Parameters.speedOfSound;
GPUParams.txFocus =UserSet.TXFocus;

GPUParams.pixelMapX = linspace(min(Trans.ElementPos(:)),max(Trans.ElementPos(:)),PData.Size(2))*1e-3;
keyboard
wavelength=GPUParams.c/Trans.frequency*1e-6/8;
GPUParams.pixelMapZ = linspace(0,PData.Size(1)*wavelength,PData.Size(1));
GPUParams.delay=(2*Trans.lensCorrection*1e-3/GPUParams.c)*ones(UserSet.na,1);
GPUParams.IQData= zeros(PData.Size(1),PData.Size(2),UserSet.numAcq,'single');
GPUParams.IQData= complex(GPUParams.IQData,0);
GPUParams.SVDEn = UserSet.SVDEn; 

% keyboard
cuFBF_int(1,GPUParams.rf_data,GPUParams.deg_tx, GPUParams.actAperture,GPUParams.filterCoef, GPUParams.pixelMapX,...
            GPUParams.pixelMapZ,GPUParams.delay,GPUParams.fs, GPUParams.ftx,GPUParams.c, GPUParams.txFocus,...
            GPUParams.reconMode,GPUParams.gpuID);
assignin('base','GPUParams',GPUParams);
return
%EF#1

%EF#2
img =GPU_beamforming(RcvData)
% keyboard;

GPUParams=evalin('base','GPUParams');
temp=permute(reshape(RcvData,GPUParams.numSample,GPUParams.na,GPUParams.numAcq,[]),[1,4,2,3]);
GPUParams.rf_data=temp(:,GPUParams.channelInd,:,GPUParams.reconFrames);
tic
img = cuFBF_int(0,GPUParams.rf_data);
toc;
GPUParams.IQData=img;
GPUParams.rf_data=temp(:,GPUParams.channelInd,:,GPUParams.reconFrames);

if GPUParams.SVDEn==1
    %Clutter Filtering
    sizeIQ=size(GPUParams.IQData);
    ns=sizeIQ(1)*sizeIQ(2);
    permInd =randperm(ns);
    
    for i=10:-1:1
        if mod(ns,i)==0
            nsplit = i;
            break;
        end
    end
    
    data=reshape(GPUParams.IQData,[],sizeIQ(3));
    permInd =reshape(permInd,[],nsplit);
    
    %%
    dataFilter = zeros(ns,sizeIQ(3),'like',data);
    for i=1:nsplit
%         tic;
%         disp(['SVD processing Part : ',num2str(i),'\',num2str(nsplit)])
        tempData=data(permInd(:,i),:);
        
        [U,S,V]=svd(tempData,'econ');
        
        %Test for distance
        logEnergy=20*log10(diag(S)/max(diag(S)));
        logEnergy1=logEnergy-min(logEnergy(:));
        logEnergy1=logEnergy1/max(logEnergy1);
        xlog=((1:length(logEnergy1))/length(logEnergy1))';
        dist=sqrt(logEnergy1.^2+xlog.^2);
        [~,index1]=min(dist);
        dataFilter(permInd(:,i),:)=U*S(:,index1:end)*V(:,index1:end)';
    end
    
    img=double(reshape(mean(dataFilter(:,1:end).*conj(dataFilter(:,1:end)),2),sizeIQ(1),sizeIQ(2)));
else
    %Adaptive filtering
%     img=double(abs(GPUParams.IQData(:,:,1)));
    img=double(mean(abs(GPUParams.IQData),3));
end
return
%EF#2

