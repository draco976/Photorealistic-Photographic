function F = quad(wname) 
%Quad wavelet filters. 
%   F = Quad(W) returns the scaling filter associated with the qmf
%   wavelet specified by the character array, 'qdN'. 
%   Possible values for N are 1. 
%
%   We have taken this from MATLAB Documentation
%
%   Copyright 2019 The MathWorks, Inc.
 
TFlem = startsWith(wname,'qd');
if ~TFlem
    error('Wavelet short name is lem followed by filter number');
end
fnum = regexp(wname,'(\d+)','match','Once');

if isempty(fnum) 
    error('Specify a filter number as 1'); 
end 

if ~isempty(fnum) 
    num = str2double(fnum);
end

tffilt = ismember(num,[1]);

if ~tffilt
    error('Filter number must be 1');
end

 
switch num 
    case 1 
% F = [... 
%         0.02807382,0.02807382;
%         -0.060944743,0.060944743;
%         -0.073386624,-0.073386624;
%          0.41472545,-0.41472545;
%          0.7973934,0.7973934;
%          0.41472545,-0.41472545;
%         -0.073386624,-0.073386624;
%         -0.060944743,0.060944743;
%          0.02807382,0.02807382 ... 
% ];

    F= QMFDesign(8,0.3,1);
 
end

%wavemngr('add','QuadratureMirrorFilter','qd',1,'1','quad')