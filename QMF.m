[N,fpts,mag,wt]=firpmord([0.395 0.605],[1 0],[0.0005 0.0005]);
[q,err]=firpm(N,fpts,mag,wt);
q(18)=q(18)+err;
h0=firminphase(q);
g0=fliplr(h0);
k=0:17;
h1=((-1).^k).*g0;
g1=-((-1).^k).*h0;
[H0,w0] = freqz(h0,1,256);
[H1,w1] = freqz(h1,1,256);
[G0,w2] = freqz(g0,1,256);
[G1,w3] = freqz(g1,1,256);
plot(w0/pi,abs(H0),'g',w1/pi,abs(H1),'r'),title('Analysis Filters');xlabel('\omega/\pi'); ylabel('Magnitude');grid;figure;plot(w2/pi,abs(G0),'b',w3/pi,abs(G1),'y'),title('Synthesis Filters');grid;xlabel('\omega/\pi'); ylabel('Magnitude');

m = input('Input Length = ');
x=randn(1,m);
v0=conv(h0,x);
u0=downsample(v0,2);
d0=upsample(u0,2);
y0=conv(d0,g0);

v1=conv(h1,x);
u1=downsample(v1,2);
d1=upsample(u1,2);
y1=conv(d1,g1);

y=y0+y1;

[X,w4] = freqz(x,1,256);
[U0,w5] = freqz(u0,1,256);
[U1,w6] = freqz(u1,1,256);
[Y,w7] = freqz(y,1,256);

subplot(2,2,1);
plot(w4/pi,abs(X),'g'),title('Gaussian Noise-Input');xlabel('\omega/\pi'); ylabel('Magnitude');grid;
subplot(2,2,2);
plot(w7/pi,abs(Y),'r'),title('Gaussian Noise-Output');xlabel('\omega/\pi'); ylabel('Magnitude');grid;
subplot(2,2,3);
plot(w5/pi,abs(U0),'b'),title('Gaussian Noise-Intermediate Signal U0');xlabel('\omega/\pi'); ylabel('Magnitude');grid;
subplot(2,2,4);
plot(w6/pi,abs(U1),'y'),title('Gaussian Noise-Intermediate Signal U1');xlabel('\omega/\pi'); ylabel('Magnitude');grid;
