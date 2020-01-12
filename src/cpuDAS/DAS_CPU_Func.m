% CPU beamforming wrapper developed by Chee Hau Leow
% Property of ULIS group, Imperial College London
% Requests for permission should be addressed to c.leow12@ic.ac.uk/mengxing.tang@ic.ac.uk

function [BF_img] = DAS_CPU_Func1(RF_data, degX_tx,degY_tx,ActElePosX,ActElePosY,...
    ActElePosZ,filterCoef, pixelMapX, pixelMapY,pixelMapZ,t0,fs,ftx,c,focus,senseCutoff,ReconMode)
%  This function read the rf data from the location and beamform the data
% Input: 
%   1)rf_data: Input rf, (sample, channel,frame)/ (sample,channel,angle,frame)
%   2) degX_tx : A vector containing the transmit steering angle in x direction (in radian)
%   3) degY_tx : A vector containing the transmit steering angle in y direction (in radian)
%   4) actApertureX:vector containing the aperture position in x direction(m)
%   5) actApertureY:vector containing the aperture position in y direction(m)
%   6) actApertureZ:A vector containing the aperture position in z direction(m)
%   7) FilterCoef: A vector containing FIR filter coefficient 
%   8) pixelMapX: A vector containing the lateral image position (m)
%   9)pixelMapY: A vector containing the elevation image position (m)
%   10) pixelMapZ:A vector containing the axial image position (m)
%   11) delay: delay that account for transmit delay, receive delay and lens correction [s]
%   12) fs: sampling frequency
%   13) ftx: transmit frequency (for apodization function calculation)
%   14) c: speed of sound
%   15) focus: transmit focal point: 0 for plane wave, -ve for diverging wave 
%   16) senseCutoff: sensitivity cutoff (0-1)
%   17) ReconMode: 0 - LRI, 1 - HRI
% Output:
% BF_img: Reconstructed image
%

[numSample,numChannel,na,numFrame]=size(RF_data);
if (na~= length(degX_tx))
    numFrame = na ;
    na =length(degX_tx);
end
RF_data =reshape(RF_data,numSample,numChannel,na,numFrame);
RF_data= double(RF_data);
imageX=length(pixelMapX);
imageY=length(pixelMapY);
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
        if (focus ==0)          %PW reconstruction
            if (ActElePosY(1)==ActElePosY(2) && ActElePosZ(1)==ActElePosZ(2)) %2D acq
                LRI(:,:,j,i)=DAS_PW(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,degX_tx(j),fs,ftx,c,t0(j),senseCutoff);
            else %3D acq
             
            end
        elseif (focus <0)       %Diverging wave reconstruction
            if all(ActElePosZ)  %Curve array
                LRI(:,:,j,i)=DAS_Diverging_curve(RF_data(:,:,j,i),ActElePosX,ActElePosZ,pixelMapX,pixelMapZ,degX_tx(j),fs,ftx,c,t0(j),focus,senseCutoff);
            else  %Phase array               
                LRI(:,:,j,i)=DAS_Diverging_flat(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,degX_tx(j),fs,ftx,c,t0(j),focus,senseCutoff);
            end
            
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

function img=DAS_Diverging_curve(RF_data,xRF,zRF,xIm,zIm,angle,fs,ftx,c,t0,focus,senseCutoff)

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

minPos = abs(xRF(1)*sin(angle));
% minPos = abs(focus*asin(xRF(1)/focus)*sin(angle));

t_tx= sqrt((fX-xgIm).^2 + (fZ-zgIm).^2)+ fZ + minPos;
for i=1:nl
    %     disp(['line:', num2str(i),'/',num2str(nl)]);
    t_rx = sqrt((zgIm-zRF(i)).^2 + (xgIm-xRF(i)).^2);
    t_bf = (t_tx + t_rx)/c;
    
    theta=(xgIm-xRF(i))./(zgIm-zRF(i))- angle;
    theta=atan(theta);
    Rx_apo=(sin(k*sin(theta)))./(k*sin(theta)).*cos(theta);  %Calculate the element sensitivity
    Rx_apo(isnan(Rx_apo))=1;
    Rx_apo(Rx_apo<senseCutoff)=0;
    
    img = img + Rx_apo.*interp1(ts, single(RF_data(:,i)), t_bf, 'linear', 0);
end
end


function img=DAS_Diverging_flat(RF_data,xRF,xIm,zIm,angle,fs,ftx,c,t0,focus,senseCutoff)

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

t_tx= sqrt((fX-xgIm).^2 + (fZ-zgIm).^2)+ fZ;
for i=1:nl
    %     disp(['line:', num2str(i),'/',num2str(nl)]);
    t_rx = sqrt(zgIm.^2 + (xgIm-xRF(i)).^2);
    t_bf = (t_tx + t_rx)/c;
    
    theta=(xgIm-xRF(i))./zgIm- angle;
    theta=atan(theta);
    Rx_apo=(sin(k*sin(theta)))./(k*sin(theta)).*cos(theta);  %Calculate the element sensitivity
    Rx_apo(isnan(Rx_apo))=1;
    Rx_apo(Rx_apo<senseCutoff)=0;
    
    img = img + Rx_apo.*interp1(ts, single(RF_data(:,i)), t_bf, 'linear', 0);
end
end


function img=DAS_PW(RF_data,xRF,xIm,zIm,angle,fs,ftx,c,t0,senseCutoff)

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
    Rx_apo(Rx_apo<senseCutoff)=0;
    
    img = img + Rx_apo.*interp1(ts, single(RF_data(:,i)), t_bf, 'linear', 0);
end
end

