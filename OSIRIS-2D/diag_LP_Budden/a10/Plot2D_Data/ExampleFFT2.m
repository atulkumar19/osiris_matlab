
%{ 

%%%%%% 2d fft starts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

Here is a quick example. I assume a random 2D image where the horizontal 
axis (columns) represents data collected in space domain and the vertical 
axis (rows) represents data collected in time domain. The two domains 
could easily be collected in the same domain. I also chose the number of 
samples collected in either domain to be different 
(i.e., a rectangular image), but the number of samples could easily be 
the same for both dimensions.

%}

clc
close all
clear all


Nx = 64; % Number of samples collected along first dimension
Nt = 32; % Number of samples collected along second dimension
dx = 1;  % Distance increment (i.e., Spacing between each column)
dt = .1; % Time increment (i.e., Spacing between each row)

x = 0 : dx : (Nx-1)*dx; % distance
t = 0 : dt : (Nt-1)*dt; % time

data_spacetime = randn(Nt,Nx); % random 2D matrix

Nyq_k = 1/(2*dx); % Nyquist of data in first dimension
Nyq_f = 1/(2*dt); % Nyquist of data in second dimension
dk = 1/(Nx*dx);   % Wavenumber increment
df = 1/(Nt*dt);   % Frequency increment

k = -Nyq_k : dk : Nyq_k-dk; % wavenumber
f = -Nyq_f : df : Nyq_f-df; % frequency

data_wavenumberfrequency = zeros(size(data_spacetime)); % initialize data

% Compute 2D Discrete Fourier Transform
for i1 = 1 : Nx
for j1 = 1 : Nt
for i2 = 1 : Nx
 for j2 = 1 : Nt
  data_wavenumberfrequency(j1,i1) = data_wavenumberfrequency(j1,i1) + ...
   data_spacetime(j2,i2)*exp(-1i*(2*pi)*(k(i1)*x(i2)+f(j1)*t(j2)))*dx*dt;
 end
end
end
end

figure;
subplot(3,1,1);
imagesc(k,f,abs(data_wavenumberfrequency));
colorbar; v = caxis;
title('2D DFT');

fft2result = fftshift(fft2(data_spacetime))*dx*dt;

subplot(3,1,2);
imagesc(k,f,abs(fft2result));
colorbar; caxis(v);
title('FFT2');

subplot(3,1,3);
imagesc(k,f,abs(fft2result-data_wavenumberfrequency));
colorbar; caxis(v);
title('Difference');

data_time = randn(Nt,1);
data_frequency = zeros(size(data_time));

%%%%%% 2d fft ends   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%}


%%%% 1D-FFT starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all


Lx=2*pi;
Ncx=200;
dx=Lx/Ncx;
xc=dx*(0:(Ncx-1))' + 0.5*dx;

yc=5*sin(2*pi*(10/Lx)*xc)+10*cos(2*pi*(5/Lx)*xc);

minKx=1/Lx;
maxKx=1/dx; % this corresponds to sampling frequency (samples per unit length)
dkx=maxKx/Ncx;

kk=dkx*(0:(Ncx-1))';
kkplot=kk(kk<(maxKx/2));

FFT_yc=fft(yc,Ncx);
powFFT=abs(FFT_yc(1:round(Ncx/2)));

figure(1)
plot(kkplot,powFFT);


figure(2)
plot(xc,yc);

%%%% 1D-FFT ends %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%