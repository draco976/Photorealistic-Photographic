function tQMFDesign (LogFile)

addpath (fullfile (pwd, '..'));

if (nargin >= 1)
  RefFile = 'tQMFDesign.ref';
  diary (LogFile);
end


N = 32;
fsb = 0.3;

% Turn off warnings for this test
warning('off', 'MATLAB:normest:notconverge');

alpha = 1;
h = QMFDesign (N, fsb, alpha);

fprintf ('Filter coefficients:\n');
h'
QMFFilterError (h)

fprintf ('\nCheck ripple computation (all should be equal)\n')
CheckRipple (h);

% Frequency response
freqz(h);

% Reference filter from Jain & Crochiere
% Results in a filter with the following characteristics
%  Stopband edge attenuation: 33 dB
%  First lobe attenuation: 44 dB
%  Passband ripple: 0.015 dB
hh = [0.46513280;
      0.13063700;
     -0.99656700E-1;
     -0.41773659E-1;
      0.53938050E-1;
      0.16805820E-1;
     -0.33077250E-1;
     -0.58240110E-2;
      0.20216010E-1;
      0.71798260E-3;
     -0.11586330E-1;
      0.12928400E-2;
      0.58649780E-2;
     -0.16349580E-2;
     -0.23388170E-2;
      0.12488120E-2];
href = [flipud(hh); hh];

fprintf ('\nReference coefficients (from Jain/Crochiere paper):\n');
href'
QMFFilterError (href)

% ----------
if (nargin >= 1)
  diary off;
  disp (' ');
  eval (['! FC ' LogFile ' ' RefFile]);
  delete (LogFile);
end

return

%---------
function CheckRipple (h)

% Check ripple computation using different approaches

N = length (h);

% Check b[n]
Neg = diag((-1).^(0:N-1));
hL = h;
gL = 2 * h;
hH = Neg * h;
gH = -2 * Neg * h;

% Direct convolution
b = 0.5 * (conv(gL, hL) + conv(gH, hH));
Sb2M = b' * b - b(N)^2;
fprintf ('Sum of squares (n ~= N-1): %g\n', Sb2M);

% Using convolution matrix
GcL = GConv (gL);
GcH = GConv (gH);
bc = 0.5 * (GcL * hL + GcH * hH);
Sbc2M = bc' * bc - bc(N)^2;
fprintf ('Sum of squares, conv. matrix (n ~= N-1): %g\n', Sbc2M);

% Alternate calculation, omitting rows in the convolution matrix
I = [(1:2:N-3) (N+1:2:2*N-3)] + 1;
GcL = GConv (gL, I);
GcH = GConv (gH, I);
bp = 0.5 * (GcL * hL + GcH * hH);
Sbp2 = bp' * bp;
fprintf ('Sum of squares (n ~= 2k & n ~= N-1): %g\n', Sbp2);

% Calculate using matrix formulation
ALL = GcL' * GcL;
ALH = GcL' * GcH;
AHL = GcH' * GcL;
AHH = GcH' * GcH;
A = 0.25 * (ALL + (ALH * Neg) + (Neg * AHL) + (Neg * AHH * Neg));
hAh = h' * A * h;
fprintf ('Sum of squares (n ~= 2k & n ~= N-1): %g\n', hAh);

return

%----------
function Gc = GConv (g, I)
% This function fills in a convolution matrix. The matrix includes a
% subset of the rows of the convolution matrix formed from the
% filter response g. The full matrix is a 2N-1 by N Toeplitz matrix
% of the form
%  G = [g[0]  0       0      ... 0      0     ]
%      [g[1]  g[0]    0      ... 0      0     ]
%      [ ...                 ... 0      0     ]
%      [g[N-2] g[N-3] g[N-4] ... g[0]   0     ]
%      [g[N-1] g[N-2] g[N-3] ... g[1]   g[0]  ]
%      [0      g[N-1] g[N-2] ... g[2]   g[1]  ]
%      [ ...                 ...              ]
%      [ 0     0      0      ... g[N-3] g[N-2]]
%      [ 0     0      0      ... g[N-2] g[N-1]]
% The index array I is used to select rows (numbered from 1 to 2N-1) from the
% full array.
N = length (g);

% Fill in the convolution matrix from the (column) vector of
% filter coefficients
% Create the first column and first row of the matrix
gcol = zeros (2*N-1, 1);
gcol(1:N) = g;
grow = zeros (1, N);
grow(1) = g(1);

% Fill in the Toeplitz matrix
GF = toeplitz (gcol, grow);

if (nargin > 1)
  Gc = GF(I,:);
else
  Gc = GF;
end

return
