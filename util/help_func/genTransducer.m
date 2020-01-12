function Trans = genTransducer(code)

%% TO ADD NEW TRANSDCUER
% To add new transducer, you can use the template below. Just key in a new
% transducer code name and modify Num_Element, Ele_pitch and that will do.
%
%	case 'L14-5'
% 		transducer.name = 'L14-5';
% 		transducer.Num_Element = 128;
% 		transducer.Ele_pitch = 3.048e-4;
% 		transducer.Ele_width = 2.65e-4;
% 		transducer.Ele_kerf = 3.5e-5;
% 		transducer.Array = [linspace(-63.5, 63.5, transducer.Num_Element)*transducer.Ele_pitch; zeros(1,transducer.Num_Element)];
% 		transducer.Array = transducer.Array';


switch code
    case 'L22-14v'
        Trans.name = 'L22-14';
        Trans.frequency = 15.625e6;
        Trans.numElement = 128;
        Trans.pitch = 0.100e-3;
        Trans.width = 0.098e-3;
        Trans.kerf = 0.002e-3;
        Trans.position = [linspace(-63.5, 63.5, Trans.numElement)*Trans.pitch; zeros(1,Trans.numElement); zeros(1,Trans.numElement)]';
        Trans.lensCorrection=1.54 * (0.05/2145 + 0.48/1147);  %in m units
        Trans.bw = [14 22]*1e6;
        Trans.connector = (1:Trans.numElement)';
        
    case 'L7-4'
        Trans.name = 'L7-4';
        Trans.frequency = 7.813e6;
        Trans.numElement = 128;
        Trans.pitch = 0.298e-3;
        Trans.width = 0.25e-3;
        Trans.kerf = 0.048e-3;
        Trans.position = [linspace(-63.5, 63.5, Trans.numElement)*Trans.pitch; zeros(1,Trans.numElement); zeros(1,Trans.numElement)]';
        Trans.lensCorrection=0.887e-3;  %in m units
        Trans.bw = [4 7]*1e6;
        Trans.connector = (1:Trans.numElement)';
        
    case 'L12-3v'
        Trans.name = 'L12-3v';
        Trans.frequency = 7.813e6;
        Trans.numElement = 192;
        Trans.pitch = 0.20e-3;
        Trans.width = 0.17e-3;
        Trans.kerf = 0.03e-3;
        Trans.position = [linspace(-95.5, 95.5, Trans.numElement)*Trans.pitch; zeros(1,Trans.numElement); zeros(1,Trans.numElement)]';
        Trans.lensCorrection=1.183e-3;
        Trans.bw = [3 12]*1e6;
        Trans.connector = [1:128 1:64]';
        
    case 'L12-5 38mm'
        Trans.name = 'L12-5 38mm';
        Trans.numElement = 192;
        Trans.pitch = 0.1979e-3;
        Trans.width = 0.1729e-3;
        Trans.kerf = 0.025e-3;
        Trans.position = [linspace(-95.5, 95.5, Trans.numElement)*Trans.pitch; zeros(1,Trans.numElement); zeros(1,Trans.numElement)]';
        Trans.lensCorrection=2.365e-3;
        Trans.bw = [5 12]*1e6;
        Trans.connector = [1:128 1:64]';
        
    case 'L12-5 50mm'
        Trans.name = 'L12-5 50mm';
        Trans.frequency = 7.813e6;
        Trans.numElement = 256;
        Trans.pitch = .1953e-3;
        Trans.width = .1703e-3;
        Trans.kerf = 0.025e-3;
        Trans.position = [linspace(-127.5, 127.5, Trans.numElement)*Trans.pitch; zeros(1,Trans.numElement);zeros(1,Trans.numElement)]';
        Trans.lensCorrection=2.365e-3;
        Trans.bw = [5 12]*1e6;
        Trans.connector = [1:128 1:128]';
        
    case 'L11-4v'
        Trans.name = 'L11-4v';
        Trans.frequency = 6.25e6;
        Trans.numElement = 128;
        Trans.pitch = 3e-4;
        Trans.width = 2.7e-4;
        % 		transducer.kerf = 0.048e-3;
        Trans.position = [linspace(-63.5, 63.5, Trans.numElement)*Trans.pitch; zeros(1,Trans.numElement);zeros(1,Trans.numElement)]';
        Trans.lensCorrection = 1.4785e-3; % in mm units; was 5 wavelengths
        Trans.bw = [4 11]*1e6;
        Trans.connector = (1:Trans.numElement)';
        
    case 'P4-1'
        Trans.name = 'P4-1';
        Trans.frequency = 2.5e6;
        Trans.numElement = 96;
        Trans.kerf= 50e-6;
        Trans.pitch = 0.2950e-3;
        % 		transducer.kerf = 0.048e-3;
        Trans.position = [linspace(-47.5, 47.5, Trans.numElement)*Trans.pitch;zeros(1,Trans.numElement);zeros(1,Trans.numElement)]' ;
        Trans.connector = [023 022 021 041 024 042 046 043 045 044 047 018 017 048 013 020 ...
            019 014 015 016 049 050 054 051 053 052 009 055 056 011 012 005 ...
            006 007 008 010 004 003 002 001 040 039 038 037 033 034 035 036 ...
            093 094 095 096 092 091 090 089 128 127 126 125 119 121 122 123 ...
            124 117 118 073 074 120 077 076 078 075 079 080 113 114 115 110 ...
            109 116 081 112 111 082 085 084 086 083 087 105 088 108 107 106];
        Trans.lensCorrection = 2.464e-3; % in mm units; was 5 wavelengths
        Trans.bw = [1 4]*1e6;
        
    case 'C5-2v'
        Trans.name = 'C5-2v';
        Trans.frequency = 3.57e6;
        Trans.numElement = 128;
        Trans.radius = 49.57e-3;
        Trans.pitch = 0.508e-3;
        Trans.kerf = 0.048e03;
        Trans.scanangle = 128*Trans.pitch/Trans.radius;
        Trans.elevationFocus = 60e-3;
        deltatheta = Trans.scanangle/128;
        firstangle = -(Trans.scanangle/2) + 0.5*deltatheta;
        Angle = firstangle:deltatheta:-firstangle;
        Trans.position(:,1) = Trans.radius*sin(Angle);
        Trans.position(:,2) = 0;
        Trans.position(:,3) = Trans.radius*cos(Angle)-Trans.radius;
        Trans.lensCorrection = 1.035e-3; % in mm units; was 5 wavelengths
        Trans.bw = [2.15 4.99]*1e6;
        Trans.connector = (1:Trans.numElement)';
        
    case 'C5-2c'
        Trans.name = 'C5-2c';
        Trans.frequency = 3.13e6;
        Trans.numElement = 128;
        Trans.radius = 50e-3;
        Trans.pitch = 0.412e-3;
        Trans.kerf = 0.048e03;
        Trans.scanangle = 128*Trans.pitch/Trans.radius;
        Trans.elevationFocus = 60e-3;
        deltatheta = Trans.scanangle/128;
        firstangle = -(Trans.scanangle/2) + 0.5*deltatheta;
        Angle = firstangle:deltatheta:-firstangle;
        Trans.position(:,1) = Trans.radius*sin(Angle);
        Trans.position(:,2) = 0;
        Trans.position(:,3) = Trans.radius*cos(Angle)-Trans.radius;
        Trans.lensCorrection = 1.035e-3; % in mm units; was 5 wavelengths
        Trans.bw = [2.2 5.5]*1e6;
        Trans.connector = (1:Trans.numElement)';
        
    case 'LPICMUS'
        Trans.name = 'LPICMUS';
        Trans.frequency = 5.208e6;
        Trans.numElement = 128;
        Trans.width= 0.27e-3;
        Trans.pitch = 0.3e-3;
        Trans.position = [linspace(-63.5, 63.5, Trans.numElement)*Trans.pitch; zeros(1,Trans.numElement); zeros(1,Trans.numElement)]';
        Trans.lensCorrection=0;
        Trans.bw = Trans.frequency * [1-0.67/2 1+0.67/2];
        Trans.connector = (1:Trans.numElement)';
end
