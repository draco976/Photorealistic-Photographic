quad = [ 0.02807382,
        -0.060944743,
        -0.073386624,
         0.41472545,
         0.7973934,
         0.41472545,
        -0.073386624,
        -0.060944743,
         0.02807382];

familyName      = 'quadraturemirrorfilter';
familyShortName = 'quad';
familyWaveType  = 1;
familyNums      = '';
fileWaveName    = 'quad.mat';

wavemngr('add',familyName,familyShortName,familyWaveType,familyNums,fileWaveName);
wavemngr('read')