%% Update -23/07/2019
% Newer version with improve structure for computing performance

%% Location of the data
clear all;
addpath(genpath('.'));

% direc='C:\Users\chl912\Documents\CUDA\testData\';
direc='./dataset/';
savedir=[direc,'result/'];mkdir(savedir);
file = dir([direc, '*.mat']);

for fileIndex=1%1:length(file)
%     close all 
    location=[direc file(fileIndex).name];   
    load(location,'UserSet','Trans');
    videoName=regexp(file(fileIndex).name,'.mat','split');
    videoName=['gpuDAS_',videoName{1}];
    vid1=VideoWriter([fullfile(savedir,videoName),'.avi']);
    vid1.FrameRate = 100;
    vid1.Quality = 50;

    %% Define recon parameter
    Recon.startFrame = 1;
    Recon.endFrame = 2;%%UserSet.numAcq*UserSet.NumFrame;
%     Recon.NumFiringPerFrame = UserSet.na;
%     Recon.skipFrame = 1; % 1 = no skip, 2 = skip 1 frame, etc
%     Recon.sumPulse=0;       %Pulse summing 0:Bmode, 1:PI all-no summing 2: PI sum-sum2pulse
%     Trans.name=UserSet.transducer;
    Recon.type  = 'LRI';		% HRI, LRI, or CRI
    [UserSet,Trans,ImagParam] = genParams(UserSet,Trans,Recon);
    if strcmp(Trans.name,'LPICMUS')
        ImagParam.delay = -abs(Trans.position(1)*sin(ImagParam.degX_tx))/ImagParam.c;
    else
        ImagParam.delay= 2*Trans.lensCorrection/ImagParam.c*ones(UserSet.na,1);
    end
    ImagParam.senseCutoff = 0.5;
    
    %%  Define filter parameters
    FiltParam.type='FIRBand';    %all,FIRHigh,FIRLow,FIRBand
    FiltParam.Fpass1=Trans.bw(1);
    FiltParam.Fpass2=Trans.bw(2);
      
    %% Define pixel map
    pixelMap.type=Trans.name(1);
    switch pixelMap.type
        case 'L'
            pixelMap.focus=0;
            pixelMap.dz=min(ImagParam.c/ImagParam.fs/2, 50e-5);
            pixelMap.dy=0;
            pixelMap.dx=Trans.pitch;
            pixelMap.upperLeft = [Trans.position(Trans.startAperture,1)*1,0, 5/1000];
            %             pixelMap.bottomRight = [Trans.position(Trans.startAperture+127,1)*1+pixelMap.dx/2 UserSet.sampleDepthMM/1000];
            pixelMap.bottomRight = [Trans.position(Trans.startAperture+127,1)*1,0, 50/1000];
        case 'C'  %Curvilinear
            aperture = Trans.numElement*Trans.pitch; % aperture based on X elements
            pixelMap.focus = -Trans.radius;
            pixelMap.dz=1540/(Trans.frequency*4*2);
            pixelMap.dy=0;
            pixelMap.dx=Trans.pitch;
            pixelMap.upperLeft = [-(UserSet.sampleDepthMM +Trans.radius)*sin(Trans.scanangle/2)*2, 0, 0*1e3]/1000;
            pixelMap.bottomRight = [(UserSet.sampleDepthMM +Trans.radius)*sin(Trans.scanangle/2)*2, 0, UserSet.sampleDepthMM+2]/1000;
        case 'P'
            aperture = Trans.numElement*(Trans.pitch); % aperture based on X elements
            pixelMap.focus = -(aperture/2)/tan(-UserSet.scanRangle);
            pixelMap.dz=1540/(Trans.frequency*4*2);
            pixelMap.dy=0;
            pixelMap.dx=Trans.pitch;
            pixelMap.upperLeft = [-UserSet.sampleDepthMM, 0, 0*1e3]/1000;
            pixelMap.bottomRight = [UserSet.sampleDepthMM, 0, UserSet.sampleDepthMM]/1000;
    end
    
    %% GPU Beamforming
    tic;
    data = DAS_CPU_beamforming(location, Trans,ImagParam, FiltParam,pixelMap,'on');
    toc
%     return
  %% Display and save Images
    clear F
    open(vid1);
    figure
    z=(0:size(data,1)-1)*pixelMap.dz*1e3+pixelMap.upperLeft(2)*1e3;
    %     x=((0:size(data,2))*pixSize+pixelMap.UpperLeft(1))*1e3;
    x=((0:size(data,2)-1)*pixelMap.dx+pixelMap.upperLeft(1))*1e3;
    
    if strcmp(ImagParam.returnType,'LRI')
        data1=squeeze(sum(data,3));
    else
        data1=data;
    end
    normMax=mean(max(max(abs(data1(:)))));
    
    for i=1:size(data1,3)
        imagesc(x,z,20*log10(abs(data1(:,:,i))/normMax),[-50,0]);
        %         imagesc(x,z,20*log10(abs(sum(data(:,:,:,i),3)/normMax)),[-30,0]);
        title(['Frame=',num2str(i)]);
        xlabel('mm');
        ylabel('mm');
        axis image;
        %         xlim([-50 50]);
        colormap(gray);
        set(gca,'FontSize',14);
        %         pause();
        F(i)=getframe(gcf);
        drawnow
        writeVideo(vid1,F(i));
    end  
    %
    close(vid1);
       
    %% save IQData
    
%     data=IQData{1};
%     clear IQData
    S=whos('data');
    splitPart=ceil(S.bytes/2^32);
    
    switch ImagParam.returnType
         case 'LRI'	
             n=4;
             if splitPart~=1
                 for i=2:splitPart
                     part=ceil((i-1)*size(data,n)/splitPart+1):ceil(i*size(data,n)/splitPart);
                     eval(['data',num2str(i),'=data(:,:,:,part);']);
                 end
                 part=1:ceil(1*size(data,n)/splitPart);
                 data=data(:,:,:,part);
             end
         case 'HRI'
             n=3;
             for i=2:splitPart
                 part=ceil((i-1)*size(data,n)/splitPart+1):ceil(i*size(data,n)/splitPart);
                 eval(['data',num2str(i),'=data(:,:,:,part);']);
             end
             part=1:ceil(1*size(data,n)/splitPart);
             data=data(:,:,part);
     end

    data=serialize(data);
    if splitPart~=1
        for i=2:splitPart
            eval(['data',num2str(i),'=serialize(data',num2str(i),');']);
        end
        
        string=['savefast([savedir,videoName,''_IQ_3_13MHz'',ImagParam.returnType],''data'',''UserSet'',''pixelMap'',''splitPart',''''];
        for i=2:splitPart
            var=['''data',num2str(i),''''];
            string=[string,',',var];
        end
        string=[string,')'];
        eval(string);
    else
        savefast([savedir,videoName,'_IQ_3_13MHz'],'data','UserSet','pixelMap');
    end
    %     
end