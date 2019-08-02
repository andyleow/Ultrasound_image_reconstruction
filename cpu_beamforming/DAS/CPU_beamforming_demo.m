% Demo function to load and reconstruct ultrafast ultrasound images 
% Written by Chee Hau Leow, ULIS group, Imperial College London
% Email : c.leow12@imperial.ac.uk

%% Load rf and settings 

addpath(genpath('C:\Users\chl912\Documents\git\cpu_beamforming\DAS\'));
load('demo_data');
focus =0;
ReconMode =0;  % 0 for LRI , 1 for HRI

%% Execute reconstruction 
BF_img = DAS_CPU_Func(rf_data(:,:,1:10), tx_angle, ActElePosX,filterCoef, pixelMapX, pixelMapZ,fs,ftx,c,focus,t0,ReconMode);

%% Display images
maxPower = max(abs(BF_img(:)));

figure,
for i=1:size(BF_img,3)
%     subplot(1,2,1)
    imagesc(pixelMapX*1e3,pixelMapZ*1e3,20*log10(abs(BF_img(:,:,i)/maxPower)),[-60 0]);
    title(['frame',num2str(i)]);
    %     image(ImgData(1:nz,:,i));
    colormap(gray);
    xlabel('mm');
    ylabel('mm');
    axis image;
    axis tight
    drawnow()
end

