% GPU beamforming wrapper developed by Chee Hau Leow
% Property of ULIS group, Imperial College London
% Requests for permission should be addressed to c.leow12@ic.ac.uk/mengxing.tang@ic.ac.uk

function [BF_img] = FBF_CPU_Func(RF_data, tx_angle, ActElePosX,filterCoef, pixelMapX, pixelMapZ,t0,fs,ftx,c,focus,ReconMode)
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

% %Filtering
% for i=1:numFrame
%      RF_data(:,:,:,i)=filter(filterCoef,1,RF_data(:,:,:,i));
% end
% t0 = t0 + 0.5*(length(filterCoef)-1)/fs;

%Delay and sum
tic;
LRI=zeros(imageZ,imageX,na,numFrame);
for i=1:numFrame
%     disp(['Frame:',num2str(i)]);
    for j=1:na
        disp(['LSR:',num2str(j),'/',num2str(na)]);      
        if (focus ==0)
            %PW reconstruction
            LRI(:,:,j,i)=FBF_PW(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,tx_angle(j),fs,ftx,c,t0(j),filterCoef);
        elseif (focus <0)
            %Diverging wave reconstruction
            LRI(:,:,j,i)=FBF_Diverging(RF_data(:,:,j,i),ActElePosX,pixelMapX,pixelMapZ,tx_angle(j),fs,ftx,c,t0(j),filterCoef,focus);
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
function img=FBF_PW(RF_data,xRF,xIm,zIm,deg_tx,fs,ftx,c,t0,FiltCoef)

[nt,nx] = size(RF_data);
ntFFT =2^(nextpow2(nt+1)); 
if length(xIm)>nx
    nxFFT=2^(nextpow2(length(xIm)+1));
else
    nxFFT =2^(nextpow2(nx+1));
end
pitch=xRF(2)-xRF(1);

%Echo Fourier space
f=(0:ntFFT/2)*fs/ntFFT; 
k=f/c;
kx=[0:nxFFT/2 -nxFFT/2+1:-1]/pitch/nxFFT;

%1d temporal fourier transfrom 
S = fft(RF_data,ntFFT,1);

%Fast time filtering
FiltCoefFFT=fft(FiltCoef,ntFFT);
S= S.*FiltCoefFFT(:);

 %Hilbert transform to get analytic signal
if 2*floor(ntFFT/2)==ntFFT %even 
    S=S(1:ntFFT/2+1,:,:);
    S(2:end-1,:,:)=2*S(2:end-1,:,:);
else %odd
    S=S(1:(ntFFT+1)/2,:,:);
    S(2:end,:,:)=2*S(2:end,:,:);   
end

%Axial delay compensation Polar-shifted to time zero 
sinA=sin(deg_tx);
tau =pitch/c*sinA*((0:nx-1)-(nx-1)*(sinA<0)) + t0+ length(FiltCoef)/2/fs;  % delay (suppose should be -ve, but sign was drop when exp(jwt) instead of exp(-jwt) was used
S=S.*exp(1i*2*pi*f'.*tau);  %kronecker tensor product of f and tau %exp(jwt) =time/phase shifted foward

%Spatial Fourier transform 
S =fft(S,nxFFT,2);

% Lateral Fourier shift pixelMap
deltaX = min(xIm(:))-min(xRF(:));
x_shift = exp(1i*2*pi*deltaX*kx);
S= S.*x_shift;

% Object fourier space
pitchImgZ = zIm(2)-zIm(1); 

% Spatial resampling base on pixelMapZ frequency 
kz=(0:ntFFT/2)/ntFFT/pitchImgZ; %Populate only positive spectrum and nyquist frequecy fs/2

%Image Fourier grid
[kzg,kxg]=ndgrid(kz,kx);

img=zeros(size(kzg));

CoefSinc = fs/(2*ftx)*pi* pitch * kxg; %Element sensitivity, originally without multiplication of 2
maskEleSense = sin(CoefSinc+eps)./(CoefSinc+eps);
maskEleSense(maskEleSense<0)=0;

%warped signal fourier space
cosA=cos(deg_tx);
sinA=sin(deg_tx);
kWarp = 0.5*(kzg.^2+kxg.^2)./(kzg*cosA+kxg*sinA+eps);

%  Create evanescent wave mask
maskEv= (kzg>=(cosA/(sinA+1)*kxg)).*(kzg>=(cosA/(sinA-1)*kxg));

for j=1:nxFFT
    img(:,j)=interp1(k,S(:,j),kWarp(:,j),'linear',0);
end

%Exclude evanescent wave
img=img.*maskEv.*maskEleSense;

%Inverse 2d FFT 
img=ifft(ifft(img,nxFFT,2),ntFFT,1);

ntshift =round(zIm(1)*2*fs/c);
ntLength =length(zIm);

img= img((1:ntLength)+ntshift,1:length(xIm));
end

function img=FBF_Diverging(RF_data,xRF,xIm,zIm,deg_tx,fs,ftx,c,t0,FiltCoef,focus)

[nt,nx] = size(RF_data);
ntFFT =2^(nextpow2(nt+1)); 
if length(xIm)>nx
    nxFFT=2^(nextpow2(length(xIm)+1));
else
    nxFFT =2^(nextpow2(nx+1));
end
pitch=xRF(2)-xRF(1);

%Echo Fourier space
f=(0:ntFFT/2)*fs/ntFFT; 
k=f/c;
kx=[0:nxFFT/2 -nxFFT/2+1:-1]/pitch/nxFFT;

%1d temporal fourier transfrom 
S = fft(RF_data,ntFFT,1);

%Fast time filtering
FiltCoefFFT=fft(FiltCoef,ntFFT);
S= S.*FiltCoefFFT(:);

 %Hilbert transform to get analytic signal
if 2*floor(ntFFT/2)==ntFFT %even 
    S=S(1:ntFFT/2+1,:,:);
    S(2:end-1,:,:)=2*S(2:end-1,:,:);
else %odd
    S=S(1:(ntFFT+1)/2,:,:);
    S(2:end,:,:)=2*S(2:end,:,:);   
end

%Axial delay compensation Polar-shifted to time zero 
sinA=sin(deg_tx);
cosA=cos(deg_tx);
fX = focus * sinA;
fZ = focus * cosA;
tau = (sqrt(fZ^2+(fX-xRF(:)').^2)+fZ)/c + t0+ length(FiltCoef)/2/fs; % delay (suppose should be -ve, but sign was drop when exp(jwt) instead of exp(-jwt) was used
S=S.*exp(1i*2*pi*f'.*tau);  %kronecker tensor product of f and tau %exp(jwt) =time/phase shifted foward

%Spatial Fourier transform 
S =fft(S,nxFFT,2);

% Lateral Fourier shift pixelMap
deltaX = min(xIm(:))-min(xRF(:));
x_shift = exp(1i*2*pi*deltaX*kx);
S= S.*x_shift;

% Object fourier space
pitchImgZ = zIm(2)-zIm(1); 

% Spatial resampling base on pixelMapZ frequency 
kz=(0:ntFFT/2)/ntFFT/pitchImgZ; %Populate only positive spectrum and nyquist frequecy fs/2

%Image Fourier grid
[kzg,kxg]=ndgrid(kz,kx);

img=zeros(size(kzg));

CoefSinc = fs/(ftx)*pi* pitch * kxg; %Element sensitivity, originally without multiplication of 2
maskEleSense = sin(CoefSinc+eps)./(CoefSinc+eps);
maskEleSense(maskEleSense<0)=0;

%warped signal fourier space
kWarp = 0.5*(kzg.^2+kxg.^2)./(kzg*cosA+kxg*sinA+eps);

%  Create evanescent wave mask
maskEv= (kzg>=(cosA/(sinA+1)*kxg)).*(kzg>=(cosA/(sinA-1)*kxg));

for j=1:nxFFT
    img(:,j)=interp1(k,S(:,j),kWarp(:,j),'linear',0);
end

%Exclude evanescent wave
img=img.*maskEv.*maskEleSense;

%Inverse 2d FFT 
img=ifft(ifft(img,nxFFT,2),ntFFT,1);

ntshift =round(zIm(1)*2*fs/c);
ntLength =length(zIm);

img= img((1:ntLength)+ntshift,1:length(xIm));
end