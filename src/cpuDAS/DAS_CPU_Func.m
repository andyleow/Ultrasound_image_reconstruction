% GPU beamforming wrapper developed by Chee Hau Leow
% Property of ULIS group, Imperial College London
% Requests for permission should be addressed to c.leow12@ic.ac.uk/mengxing.tang@ic.ac.uk

function [BF_img] = DAS_CPU_Func(RF_data, tx_angle, ActElePosX,filterCoef, pixelMapX, pixelMapZ,t0,fs,ftx,c,focus,ReconMode)
%  This function read the rf data from the location and beamform the data
% Input:
%   1)rf_data: Input rf, multiPW(sample, channel,angle,frame)/singlePW(sample, channel, frame)
%   2) tx_angle : A vector containing the transmit steering angle (in radian)
%   3) ActElePosX:A vector containing the aperture position
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
%

[numSample,numChannel,na,numFrame]=size(RF_data);
if (na~= length(tx_angle))
    numFrame = na ;
    na =length(tx_angle);
end
RF_data =reshape(RF_data,numSample,numChannel,na,numFrame);
RF_data= double(RF_data);
imageX=length(pixelMapX);
imageZ=length(pixelMapZ);

%Filtering
for i=1:numFrame
     RF_data(:,:,:,i)=filter(filterCoef,1,RF_data(:,:,:,i));
end
t0 = t0 + 0.5*(length(filterCoef)-1)/fs;

%Hilbert Transform
for i=1:numFrame
    RF_data(:,:,:,i)=hilbert(RF_data(:,:,:,i));
end

%Delay and sum
tic;
LRI=zeros(imageZ,imageX,na,numFrame);
for i=1:numFrame
%     disp(['Frame:',num2str(i)]);
    for j=1:na
        disp(['LSR:',num2str(j),'/',num2str(na)]);      
        if (focus ==0)
            %PW reconstruction
            LRI(:,:,j,i)=DAS_PW(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,tx_angle(j),fs,ftx,c,t0(j));
        elseif (focus <0)
            %Diverging wave reconstruction
            LRI(:,:,j,i)=DAS_Diverging(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,tx_angle(j),fs,ftx,c,t0(j),focus);
        end
    end
end
toc;

if ReconMode ==1
    BF_img=LRI;
else
    HRI = squeeze(sum(LRI,3));
    BF_img=HRI;
end

end

function img=DAS_Diverging(RF_data,xRF,xIm,zIm,angle,fs,ftx,c,t0,focus)
%     image=zeros(PixelMap.X,PixelMap.Z);

%Image Pixel position
%Focal point coordinate
fX = focus * sin(angle);
fZ = focus * cos(angle);

[xgIm,zgIm]=meshgrid(xIm,zIm);
img = zeros(size(xgIm));
[ns,nl] = size(RF_data);

ts = (linspace(0,ns-1,ns)/fs-t0)';

d=xRF(2)-xRF(1);
lambda=c/(ftx);
k=pi*d/lambda;

for i=1:nl
%     disp(['line:', num2str(i),'/',num2str(nl)]);
    t_tx= sqrt((fX-xgIm).^2 + (fZ-zgIm).^2)+ fZ;
    t_rx = sqrt(zgIm.^2 + (xgIm-xRF(i)).^2);
    t_bf = (t_tx + t_rx)/c;
    
    theta=(xgIm-xRF(i))./(zgIm)- angle;
    theta=atan(theta);
    Rx_apo=(sin(k*sin(theta)))./(k*sin(theta)).*cos(theta);  %Calculate the element sensitivity
    Rx_apo(isnan(Rx_apo))=1;

    img = img + Rx_apo.*interp1(ts, single(RF_data(:,i)), t_bf, 'linear', 0);
end
end


function img=DAS_PW(RF_data,xRF,xIm,zIm,angle,fs,ftx,c,t0)

%Image Pixel position
[xgIm,zgIm]=meshgrid(xIm,zIm);

[ns,nl] = size(RF_data);
img= zeros(size(xgIm));

minPos = abs(xRF(1)*sin(angle));

ts = (linspace(0,ns-1,ns)/fs-t0)';

d=xRF(2)-xRF(1);
lambda=c/(ftx);
k=pi*d/lambda;

for i=1:nl
    %         disp(['line:', num2str(i),'/',num2str(nl)]);
    t_tx= (zgIm*cos(angle) + (xgIm*sin(angle)+ minPos));
    t_rx = sqrt(zgIm.^2+(xgIm-xRF(i)).^2);
    t_bf = (t_tx + t_rx)/c;
    
    theta=(xgIm-xRF(i))./(zgIm) - angle;
    theta=atan(theta);
    Rx_apo=(sin(k*sin(theta)))./(k*sin(theta)).*cos(theta);  %Calculate the element sensitivity
    Rx_apo(isnan(Rx_apo))=1;
    img = img + Rx_apo.*interp1(ts, single(RF_data(:,i)), t_bf, 'linear', 0);
end
end



% function img=DAS_PW(RF_data,xRF,xIm,zIm,angle,fs,ftx,c,t0)
%
%     %Image Pixel position
%     [xgIm,zgIm]=meshgrid(xIm,zIm);
%
%     [ns,nl] = size(RF_data);
%     [nr2,nc2] =deal(length(zIm),length(xIm));
%
%      minPos = abs(xRF(1)*sin(angle));
%
%     zRF = (linspace(0,ns-1,ns)/fs-t0)'*c/2;
%     ts = zRF*2/c;
%
%     tempImg=zeros(nr2,nc2,nl);
%     d=xRF(2)-xRF(1);
%     lambda=c/(ftx);
%     k=pi*d/lambda;
%
%     for i=1:nl
% %         disp(['line:', num2str(i),'/',num2str(nl)]);
%         t_tx= (repmat(zIm*cos(angle),1,nc2) + repmat(xIm*sin(angle) + minPos,nr2,1))/c;
%         t_rx = sqrt(repmat(zIm.^2,1,nc2)+repmat((xIm-xRF(i)).^2,nr2,1))/c;
%         t_bf = t_tx + t_rx;
%
%         theta=(xgIm-xRF(i))./(zgIm);
%         theta=atan(theta);
%         Rx_apo=(sin(k*sin(theta)))./(k*sin(theta)).*cos(theta);  %Calculate the element sensitivity
%         Rx_apo(isnan(Rx_apo))=0;
% %         Rx_apo=interp2(xgApo,zgApo,Apo,xgIm-xApo(i),zgIm,'linear',0);
%         tempImg(:,:,i) = Rx_apo.*interp1(ts, single(RF_data(:,i)), t_bf, 'linear', 0);
%     end
%     img=sum(tempImg,3);
% end
