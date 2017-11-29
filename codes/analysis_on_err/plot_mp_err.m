%  sndemo33.m    Multipath Analysis
%
clear all
close all

path = './StaNAV_commercial_code/satnav3p04';
addpath(path)

%%

delta = mpgen(24,3600,1,33333);
figure;
plot (delta (1:60, 1),'r','linewidth',2)
title('Multipath Error received by the 8th satellite');
xlim([0,60]);
ylim([-6,8]);
ylabel('residuals in meters')
xlabel('run time in seconds')
legend('multipath error')
%%
% Try plot multipath error
load days.mat
file_path = './StaNAV_commercial_code/satnav3p04';
addpath(file_path)

mpmat=mpgen(24,3600,1,33333);
refllh = [0*pi/180 0*pi/180 0];
refxyz = llh2xyz(refllh);

loadgps
i=0;
bar1 = waitbar(0,'Generating Ranges...  ');
randn('state',74347098);
for t = 41001:1:41060
    i=i+1;
    [svmatref,svidref] = gensv(refxyz,t,0);
    [prref,adrref] = genrng(1,refxyz,svmatref,svidref,t,...
                            [1 1 0 1 1],[],mpmat);
    pr8(i) = prref(find(svidref==8));
    adr8(i) = adrref(find(svidref==8));
	 waitbar(i/60)    
end
close(bar1); 
resi = pr8 - adr8;
% resi = resi - mean(resi);
figure;
plot(resi)
title('Multipath, Thermal Noise and Iono Divergence')
ylabel('residuals in meters')
xlabel('run time in seconds')

%% 
% The following section is for solnly ploting multipath error
x = 1:60;
p = polyfit(x,resi,1);
yfit = polyval(p,x);
figure;
plot(x,resi,x,yfit)
title('Residuals and Straight-Line Fit')

figure;
resi_fit = resi - yfit;
plot(x,resi_fit)
title('Multipath and Noise with Iono Removed')
ylabel('residuals in meters')
xlabel('run time in seconds')
%%
% hoursPerDay = 24;
% coeff24hMA = ones(1, hoursPerDay)/hoursPerDay;
P=6;
coeff = ones(1, P)/P;
% avg24hTempC = filter(coeff24hMA, 1, tempC);
avg24h = filter(coeff, 1, resi_fit);
figure;
plot(resi_fit);
hold on;
plot(avg24h,'linewidth',2);
% hold on;
% plot(delta (1:600,8),'linewidth',2)
legend('multipath error with noise','medium frequency noise','True Multipath error','location','best')
ylabel('residuals in meters')
xlabel('run time in seconds')
title('Multipath error, with Noise and Iono Removed')
%%

% x = [0:0.1:2*pi];
% y = sin(x); 

%% 
% Long time constant autoregressive model
b = [0.01087642013487   0.01087642013487];
a = [1.00000000000000  -0.97824715973025];
      
x = 15*randn(60,1);
y = filter(b,a,x);

% fft: Discrete Fourier transform.
z = fft(y); 
z = fftshift(z); 
N = length(y); 
f = [-N/2:N/2-1]/N; 
figure;
plot(f,abs(z),'o-.')
grid;
title('Discrete Fourier transform long term autoregressive model')

figure;
freqz(b,a,60)
title('Freqz long term autoregressive model')

%%
% Short time constant autoregressive model
b = [0.03780475417090   0.03780475417090];
a = [1.00000000000000  -0.92439049165821];

     
x = 15*randn(60,1);
y = filter(b,a,x);

% fft: Discrete Fourier transform.
z = fft(y); 
z = fftshift(z); 
N = length(y); 
f = [-N/2:N/2-1]/N; 
figure;
plot(f,abs(z),'o-.')
grid;
title('Discrete Fourier transform short term autoregressive model')

figure;
freqz(b,a,60)
title('Freqz short term autoregressive model')
