% GPU beamforming wrapper developed by Chee Hau Leow
% Property of ULIS group, Imperial College London
% Requests for permission should be addressed to c.leow12@ic.ac.uk/mengxing.tang@ic.ac.uk 

function [BF_img] = DAS_CPU_Func(RF_data, tx_angle, ActElePosX,filterCoef, pixelMapX, pixelMapZ,fs,ftx,c,focus,t0,ReconMode)
%  This function read the rf data from the location and beamform the data
% Input: 
%   1)rf_data: Input rf, multiPW(sample, channel,angle,frame)/singlePW(sample, channel, frame)
%   2) tx_angle : A vector containing the transmit steering angle (in radian)
%   3) ActElePosX:A vector containing the aperture position
%   4) filterCoef: A vector containing FIR filter coefficient
%   5) pixelMapX: A vector containing the lateral image position
%   6) pixelMapZ:A vector containing the axial image position
%   7) fs: sampling frequency
%   8) ftx: transmit frequency (for apodization function calculation)
%   9)  c: speed of sound
%   10) focus: transmit focal point: 0 for plane wave, -ve for diverging wave 
%   11) t0: delay that account for transmit delay, receive delay and lens correction [s]
%   12)  ReconMode: 0 for LRI and 1 for HRIclear all
% Output:
% BF_img: Reconstructed image
%

[numSample,numChannel,na,numFrame]=size(RF_data);
if (na~= length(tx_angle))
    numFrame = na ;
    na =length(tx_angle);
end
RF_data =reshape(RF_data,numSample,numChannel,na,numFrame);
   
    
imageX=length(pixelMapX);
imageZ=length(pixelMapZ);

%Filtering
for i=1:numFrame
    RF_data(:,:,:,i)=filter(filterCoef,1,RF_data(:,:,:,i));
end

%Hilbert Transform
for i=1:numFrame
    RF_data(:,:,:,i)=hilbert(RF_data(:,:,:,i));
end

%Delay and sum
tic;
LRI=zeros(imageZ,imageX,na,numFrame);
for i=1:numFrame
    disp(['Frame:',num2str(i)]); 
    for j=1:na
%         disp(['LSR:',num2str(j)]);  
%             LRI(:,:,j,i)=DAS_beamforming_PW_dynamicApo(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,angle(j),Rx_apo,fs,c,t0(j));

        if (focus ==0)
            %PW reconstruction
            LRI(:,:,j,i)=DAS_beamforming_PW_dynamicApo1(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,tx_angle(j),fs,ftx,c,t0(j));
        elseif (focus <0)
            %Diverging wave reconstruction
        end
            
            
    end
end
toc;

if strcmp(ReconMode,'LRI')
    BF_img=LRI;
else
    HRI = squeeze(sum(LRI,3));
    BF_img=HRI;
end

end

function img=DAS_beamforming_PW_dynamicApo1(RF_data,xRF,xIm,zIm,angle,fs,ftx,c,t0)
%     image=zeros(PixelMap.X,PixelMap.Z);
    
    %Image Pixel position
    [xgIm,zgIm]=meshgrid(xIm,zIm);
   
    [ns,nl] = size(RF_data);
    [nr2,nc2] =deal(length(zIm),length(xIm));
    
%     z = (linspace(0,nl-1,nl)/Trans.fs + (Trans.t0))'*Trans.c0/2;
    zRF = (linspace(0,ns-1,ns)/fs-t0)'*c/2;
    ts = zRF*2/c;
     
    if angle <0
        minPos = xRF(1)*sin(angle);
    else
        minPos = -xRF(1)*sin(angle);
    end
    
    tempImg=zeros(nr2,nc2,nl);
    d=xRF(2)-xRF(1);
    lambda=c/(ftx);
    k=pi*d/lambda;
    
    for i=1:nl
%         disp(['line:', num2str(i),'/',num2str(nl)]);
        t_tx= (repmat(zIm*cos(angle),1,nc2) + repmat(xIm*sin(angle) + minPos,nr2,1))/c;
        t_rx = sqrt(repmat(zIm.^2,1,nc2)+repmat((xIm-xRF(i)).^2,nr2,1))/c;
        t_bf = t_tx + t_rx;
       
        theta=(xgIm-xRF(i))./(zgIm);
        theta=atan(theta);
        Rx_apo=(sin(k*sin(theta)))./(k*sin(theta)).*cos(theta);  %Calculate the element sensitivity
        Rx_apo(isnan(Rx_apo))=0;
%         Rx_apo=interp2(xgApo,zgApo,Apo,xgIm-xApo(i),zgIm,'linear',0);
        tempImg(:,:,i) = Rx_apo.*interp1(ts, single(RF_data(:,i)), t_bf, 'linear', 0);
    end
    img=sum(tempImg,3);
end