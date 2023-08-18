function h = QMFDesign (N, fsb, alpha)
% Design a two-channel QMF filter bank
% The four filters in the filter bank are based on a lowpass prototype H(z)
% which is designed by this procedure. The analysis filters are HL(z) and
% HH(z) and the synthesis filters are GL(z) and GH(z),
%   HL(z) = H(z);      HH(z) = H(-z)
%   GL(z) = 2 H(z);    GH(z) = -2 H(-z)
% This filter bank has an alias cancelling property, but does introduce
% a ripple into the end-to-end frequency response. The design procedure
% minimizes a weighted combination of the ripple energy and the stopband
% energy of the lowpass prototype filter. The design procedure is motivated
% by a paper by Jain and Crochiere.
%   V. K. Jain and R. E. Crochiere, "Quadrature Mirror Filter
%   Design in the Time Domain", IEEE Trans. Acoustics, Speech,
%   Signal Processing, vol. 32, no. 2, pp. 253-361, April 1984.
%
% N - Number of coefficients for the lowpass prototype. N must be even.
% fsb - Normalized stopband edge frequency for the lowpass prototype, where
%     0.25 < fsb < 0.5.
% alpha - Relative weighting between the stopband energy and the ripple
%     in the overall response. Increasing alpha will result in greater
%     stopband attenuation in the lowpass prototype.
%
% This routine returns the coefficients of the symmetric prototype filter
% H(z). The coefficients are normalized to have energy 1/2. When the
% analysis and synthesis filters are configured as shown above, the overall
% system has an impulse response with a unity coefficient at a delay of
% N-1 samples.
%
% A sample filter can be designed with QMFDesign (32, 0.3, 1);

% $Id: QMFDesign.m,v 1.13 2009/06/03 12:21:33 pkabal Exp $

if (mod (N, 2) ~= 0)
  error ('QMFDesign: No. coefficients must be even');
end
M = N/2;

% System under consideration (ignoring aliasing)
%   B(z) = 1/2 ( GL(z) HL(z) + GH(z) HH(z) )
% Each filter has N = 2*M coefficients and is symmetric.
% The component filters are based on a symmetric lowpass prototype,
%   H(z) = z^(-(N-1)) H(1/z)
% The polyphase decomposition of H(z) is
%   H(z) = H0(z^2) + z^(-1) H1(z^2)
% where because of the symmetry of H(z) (and N being even)
%   H1(z) = z^(-(M-1)) H0(1/z).
%
% The design strategy is to find H0(z), use this to form H(z) and then get
% the component filters. An iterative procedure is used.
% - Start with an initial H0(z)
% - Begin iteration
%   - Find H(z) from H0(z).
%   - Find GL(z) and GH(z) from H(z)
%   - Express the B(z) in terms of GL(z), GH(z), HL(z) and HH(z), where
%     HL(z) and HH(z) are unknown, but depend only on a new H0(z) which is
%     to be determined.
%   - Find a new H0(z) which minimizes the ripple in B(z) and the stopband
%     energy in H(z). Note that this H0(z) will not necessarily cancel
%     aliasing.
%   - Iterate
% In the limit, there will be no change in H(z). Using this value of H(z),
% the filters in the alias cancelling filter bank can be determined.

% Notes:
% - One can use a permutation of H0(z) in the design. One such permutation
%   is the first (or second) half of of the symmetric lowpass filter.
% - The stopband energy criterion also affects the passband ripple of
%   the lowpass prototype. The ripple criterion (on the overall filter
%   bank impulse response) essentially forces the transition region
%   (pi-fsb, fsb) of the prototype to have the correct symmetry such that
%   frequency response of the overall filter bank is (approximately)
%   constant.
% - The ripple error measurement omits terms known to be zero in the
%   overall filter bank response. The design procedure gives the "same"
%   results if these terms are included.
% - This implementation of the Jain & Crochiere algorithm does not
%   converge when implemented as per their paper. The error
%   (h0' (PAP + alpha PFP) h0) does not decrease monotonically.
%   The fix was a more cautious update:
%     - Previous filter is h0p
%     - New filter via eigen solution is h0n
%     - Updated filter is h0 = (1-t) h0n + t h0p

% Initialization of the lowpass prototype
h = zeros (N, 1);
h(N/2) = 1/2;
h(N/2+1) = 1/2;

% Permutation matrix
% Maps zeroth polyphase component to the prototype filter
%   h = P h0,
% where h is Nx1, P is N x M, and h0 is Mx1.
P = zeros (N, M);
P(1:2:N,:) = eye(M);
P(2:2:N,:) = flipud(eye(M));
% P = [flipud(eye(M)); eye(M)];  % h0 is half of h

% Stopband error integral matrix (Square error is h' F h)
F = StopbandE (2 * pi * fsb, N);

% Convert to a matrix operating on the polyphase component
% The stopband error is h0' (P' F P) h0.
PFP = P' * F * P;

% Main iteration loop
Niter = 30;
hp = h;
tau = 0.2;
Converge = 0;
tol = 1e-10;
for (kk = 1:Niter)

  hp = h;  % Save the current estimate

  % Form the ripple energy matrix
  A = RippleE (h);
  PAP = P' * A * P;

  % Solve the eigenproblem
  % Find the eigenvector corresponding to the smallest eigenvalue
  [h0, lambda] = svds (PAP + alpha * PFP, 1, 0);

  % Choose the sign of the coefficients to make the largest coefficient
  % positive
  % Normalize h0 to energy 1/4 (h then has energy 1/2)
  h0 = NormE (h0, 1/4);
  hn = P * h0;
  
  % Update
  h = (1 - tau) * hn + tau * hp;
  h = NormE (h, 1/2);

  % Test for convergence
  Ed = (h - hp)' * (h - hp);
  if (Ed < tol)
    Converge = 1;
    break
  end

end

if (~ Converge)
  fprintf ('QMFDesign: Maximum iterations reached');
end

return

%----------
function F = StopbandE (wsb, N)
% Stopband error integral matrix
% The frequency response of the lowpass prototype can be written as,
%   H(w) = e'(w) h , where e(w) is a vector with elements exp(-j n w).
% Then the square magnitude response is
%   |H(w)|^2 = h' e(w) e'(w) h.
% The integral of the square magnitude response from wa to pi is
%   Es(wa,pi) = h' F(wa,pi) h.

F = zeros(N, N);
for (k = 0:N-1)
  for (m = 0:N-1)
    F(k+1,m+1) = ExpInt (k-m, wsb, pi);
  end
end

return

%----------
function A = RippleE (h)
% Ripple error energy matrix

N = length (h);

% Negation matrix, creates H(-z) from H(z) or vice-versa
% Negates every second coefficient
Neg = diag((-1).^(0:N-1));

gL = 2 * h;
gH = - 2 * Neg * h;

% Form the ripple matrix
% The impulse response the system (ignoring aliasing) is
%   b[n] = 1/2 ( SUM hL[k] gL[n-k] + SUM hH[k] gH[n-k] )
%                 k                   k
% b[n] is of length 2*N-1. The convolutions can be expressed as a
% matrix-vector product,
%   b = 1/2 (GcL hL + GcH hH),
% where GcL and GcH are convolution matrices of size 2N-1 by N. Each
% row of the convolution matrices corresponds to a time index n.
%
% For the QMF filter bank, b[n] is zero for n even and b[N-1] = 1.
% The latter assumes normalization of the energy of h to 1/2. The
% optimization will be over the unconstrained coefficients. Label the
% index set of the unconstrained coefficients as I.
%   b(I) = 1/2 (GcL(I) hL + GcH(I) hH).
% The sum of squares can be expressed as
%   SUM b^2[n] = b(I)' b(I)
%   n==I
%              = 1/4 ( hL' GcL(I)' GcL(I) hL + hL' GcL(I)' GcH(I) hH
%                    + hH' GcH(I)' GcL(I) hL + hH' GcH(I)' GcH(I) hH )
% We can see four terms of the form,
%   hx' Axy hy
% where
%   Axy = Gcx(I)' Gcy(I)
% After calculating the matrices, the sum of squares is given by
%   B2 = 1/4 (hL' ALL hL) + (hL' ALH hH) + (hH' AHL hL) + (hH' AHH hH).
% HH(z) = HL(-z). Then
%   hL = h; hH = Neg h.
% With these substitutions
%  B2 = 1/4 h' (ALL + ALH Neg + Neg AHL + Neg AHH Neg) h
%     = h' A h.

% Create the correlation matrix as the product of convolution matrices.
% The rows of Gcx correspond to the index n of b[n]. In the final
% configuration, b[n] = 0, for n even, and b[N-1] = 1 for any
% prototype filter. The known coefficients are,
%   K = [0, 2, ... , N-4, N-2, N-1, N, N+2, ... , 2N-4, 2N-2]
%       [0, 0, ... , 0,   0,   1,   0, 0,   ... , 0,    0   ]
% We will then only evaluate the ripple for the remaining coefficients:
%   I = [1, 3, ... , N-5, N-3, N+1, N+3, ... , 2N-5, 2N-3]
I = [(1:2:N-3) (N+1:2:2*N-3)] + 1;
GcL = GConv (gL, I);
GcH = GConv (gH, I);

% Form the overall matrix
A = 0.25 * (GcL' * GcL + (GcL' * GcH) * Neg ...
            + Neg * (GcH' * GcL) + Neg * (GcH' * GcH) * Neg);

return

%----------
function Gc = GConv (g, I)
% This function fills in a convolution matrix. The matrix includes a subset
% of the rows of the convolution matrix formed from the filter response g.
% The full matrix is a 2N-1 by N Toeplitz matri of the form
%  G = [g[0]  0       0      ... 0      0     ]
%      [g[1]  g[0]    0      ... 0      0     ]
%      [ ...                 ... 0      0     ]
%      [g[N-2] g[N-3] g[N-4] ... g[0]   0     ]
%      [g[N-1] g[N-2] g[N-3] ... g[1]   g[0]  ]
%      [0      g[N-1] g[N-2] ... g[2]   g[1]  ]
%      [ ...                 ...              ]
%      [ 0     0      0      ... g[N-3] g[N-2]]
%      [ 0     0      0      ... g[N-2] g[N-1]]
% The index array I is used to select rows (numbered from 1 to 2N-1) from
% the full array.
N = length (g);

% Fill in the convolution matrix from the (column) vector of filter
% coefficients
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

%----------
function yI = ExpInt (k, wa, wb)
% Evaluate the integral
%  yI = 1/(2 pi) [Int(wa,wb) exp(-j w k) dw 
%              +  Int(-wb,-wa) exp(-j w k) dw]
%     = 1/(-j 2 pi k) [exp(-j wb k) - exp(-j wa k)
%                      + exp(j wa k) - exp(j wb k)]
%     = 1/(pi k) [sin(wb k) - sin(wa k)]

if (k == 0)
  yI = (wb - wa) / pi;
else
  yI = 1/(pi*k) * (sin(wb*k) - sin(wa*k));
end

return

%----------
function hN = NormE (h, E)
% Normalize h to energy E (default 1) and make the largest coefficient
% positive

if (nargin == 1)
  E = 1;
end
% Make the largest coefficient positive
[hmax, imax] = max (abs(h));
if (h(imax) < 0)
  h = -h;
end

% Normalize h to energy E
Eh = h' * h;
hN = h * sqrt(E / Eh);

return
