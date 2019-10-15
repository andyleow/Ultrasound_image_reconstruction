function Hd = computeFIR(fs,filtParam)
switch (filtParam.type)
    case 'FIRLow'
        filtParam.Fpass1=0;
        filtParam.Coef = FIRCustom(fs,filtParam);
    case 'FIRHigh'
        filtParam.Fpass2=ImagParam.fs/2;
        filtParam.Coef = FIRCustom(fs,filtParam);
    case 'FIRBand'
        filtParam.Coef = FIRCustom(fs,filtParam);
    case 'all'
        filtParam.Coef=1;
end

Fpass1 = filtParam.Fpass1;         % First Passband Frequency
Fpass2 = filtParam.Fpass2;         % Second Passband Frequency

if Fpass2>=fs/2
    Hd = fir1(51,Fpass1/(fs/2),'high');
else
    Hd = fir1(51,[Fpass1, Fpass2]/(fs/2));
end

function Hd = FIRCustom(Fs,FiltParam)

Fpass1 = FiltParam.Fpass1;         % First Passband Frequency
Fpass2 = FiltParam.Fpass2;         % Second Passband Frequency

if Fpass2>=Fs/2
    Hd = fir1(51,Fpass1/(Fs/2),'high'); 
else
    Hd = fir1(51,[Fpass1, Fpass2]/(Fs/2)); 
end