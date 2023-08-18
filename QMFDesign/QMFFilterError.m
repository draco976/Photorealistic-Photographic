function QMFFilterError (h)
% Given a prototype lowpass filter used to form a QMF filter bank, find the
% maximum ripple in the overall response of the filter bank and the minimum
% stopband attenuation of the prototype filter.
% h - Vector of filter coefficients (symmetric)

% $Id: QMFFilterError.m,v 1.4 2009/06/03 12:21:53 pkabal Exp $

h = h(:);
if (h ~= flipud(h))
  error ('QMFFilterError: Filter is not symmetric');
end

N = length (h);
if (mod(N,2) ~= 0)
  error ('QMFFilterError: No. coefficients must be even');
end

% Form the QMF component filters
Neg = diag((-1).^(0:N-1));
hL = h;
gL = 2 * h;
hH = Neg * h;
gH = -2 * Neg * h;

% Frequency samples
Nfreq = 2001;
w = linspace (0, pi, Nfreq);

% Direct calculation of the overall response b[n]
b = 0.5 * (conv(gL, hL) + conv(gH, hH));
Bw = freqz (b, 1, w);

% Overall response ripple
Delay = N - 1;
Aw = exp(1j * w * Delay) .* Bw;
if (imag(Aw) > 1e-13)
  error ('Non-negligible imaginary part');
end
Aw = real(Aw);

MaxdB = max (abs (20 * log10 (Aw)));
fprintf ('Maximum overall ripple: %g dB\n', MaxdB);

% Prototype filter stopband attenuation
Hw = freqz (h, 1, w);

% Search for the stopband edge
% We will follow the falling edge of the transition band until
% the response falls below the first sidelobe
for (i = fix(Nfreq/2):Nfreq)
  EdgeAtt = abs(Hw(i));
  StopBandAtt = max (abs (Hw(i:Nfreq)));
  if (EdgeAtt < StopBandAtt)
    break
  end
  fsb = w(i) / (2 * pi);
end

StopBandMindB = -20 * log10 (StopBandAtt);
fprintf ('Minimum stopband attenuation: %g dB\n', StopBandMindB);
fprintf ('Stopband edge: %g\n', fsb);

return
