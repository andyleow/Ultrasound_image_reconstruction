%% Update -23/07/2019
% Newer version with improve structure for computing performance

%% Location of the data
clear all;
addpath(genpath('.'));

% direc='C:\Users\chl912\Documents\CUDA\testData\';
direc='./dataset/'
% direc='/home/chleow/HAIFU_testData/';
savedir=[direc,'result/'];mkdir(savedir);
file = dir([direc, '*.mat']);

for fileIndex=8%:length(file)
    %     close all
    location=[direc file(fileIndex).name];
    load(location,'UserSet','Trans');
    videoName=regexp(file(fileIndex).name,'.mat','split');
    videoName=['gpuDAS_',videoName{1}];
    vid1=VideoWriter([fullfile(savedir,videoName),'.avi']);
    vid1.FrameRate = 100;
    vid1.Quality = 50;
    
    %% Determine the type of transducer
    %         Trans.name ='LPICMUS';
    %     Trans.name = 'P4-1';
    %         Trans.name = 'L12-3v';
    %% Define recon parameter
    Recon.startFtrame = 1;
    Recon.endFrame =1;% UserSet.numAcq*UserSet.NumFrame;
    %     Recon.NumFiringPerFrame = UserSet.na;
    %     Recon.skipFrame = 1; % 1 = no skip, 2 = skip 1 frame, etc
    %     Recon.sumPulse=0;       %Pulse summing 0:Bmode, 1:PI all-no summing 2: PI sum-sum2pulse
%     Trans.name=UserSet.tranducer;
    Recon.type  = 'HRI';		% HRI, LRI, or CRI
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
    FiltParam.tempMedFilt=0; %temporal median filter
    
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
        case 'M'
            pixelMap.focus=UserSet.focusMM*1e-3;
            pixelMap.dz=50e-6;%1540/(Trans.frequency);
            pixelMap.dx=Trans.pitch/2;
            pixelMap.dy=Trans.pitch/2;
            pixelMap.upperLeft = [min(Trans.position(:,1))*1e3-UserSet.sampleDepthMM * tand(UserSet.viewAngle),...
                min(Trans.position(:,2))*1e3-UserSet.sampleDepthMM * tand(UserSet.viewAngle),...
                0*1e3]/1000;
            pixelMap.bottomRight = [max(Trans.position(:,1))*1e3+UserSet.sampleDepthMM * tand(UserSet.viewAngle),...
                max(Trans.position(:,2))*1e3+UserSet.sampleDepthMM * tand(UserSet.viewAngle),...
                UserSet.sampleDepthMM]/1000;
    end
    
    %% GPU Beamforming
    tic;
    data = DAS_GPU_beamforming(location, Trans,ImagParam, FiltParam,pixelMap,'on');
    toc
    %     return
    %% Display and save Images
    open(vid1);
    
    figure
    z=(0:size(data,1)-1)*pixelMap.dz*1e3+pixelMap.upperLeft(2)*1e3;
    %     x=((0:size(data,2))*pixSize+pixelMap.UpperLeft(1))*1e3;
    x=((0:size(data,2)-1)*pixelMap.dx+pixelMap.upperLeft(1))*1e3;
    
    if strcmp(ImagParam.returnType,'LRI')
        if pixelMap.type=='M'
            data1=squeeze(sum(data,4));
        else
            data1=squeeze(sum(data,3));
        end
    else
        data1=data;
    end
    normMax=mean(max(max(abs(data1(:)))));
    
    if pixelMap.type=='M'
        if ndims(data1)>3  %Display centre slices for display i
            for i=1:size(data1,4)
                subplot(1,2,1)
                imagesc(x,z,20*log10(abs(squeeze(data1(:,:,round(size(data1,3)/2),i)))/normMax),[-60,0]);
                title(['SliceX, Frame=',num2str(i)]);
                xlabel('mm');
                ylabel('mm');
                axis image;
                colormap(gray);
                set(gca,'FontSize',14);
                subplot(1,2,2)
                imagesc(x,z,20*log10(abs(squeeze(data1(:,round(size(data1,2)/2),:,i)))/normMax),[-60,0]);
                title(['Slice Y, Frame=',num2str(i)]);
                xlabel('mm');
                ylabel('mm');
                axis image;
                colormap(gray);
                set(gca,'FontSize',14);
                pause(0.01)
                F(i)=getframe(gcf);
                drawnow
                writeVideo(vid1,F(i));
            end
        else
            for i=1:size(data1,3)
                imagesc(x,z,20*log10(abs(data1(:,:,i))/normMax),[-60,0]);
                title(['Slice = ',num2str(i)]);
                xlabel('mm');
                ylabel('mm');
%                 axis image;
                colormap(gray);
                set(gca,'FontSize',14);
                pause(0.01)
                F(i)=getframe(gcf);
                drawnow
                writeVideo(vid1,F(i));
            end
        end
    else
        for i=1:size(data1,3)
            imagesc(x,z,20*log10(abs(data1(:,:,i))/normMax),[-60,0]);
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
    end
    %
    close(vid1);
    
    %     continue
        return
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