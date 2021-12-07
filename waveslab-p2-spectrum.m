close all
clear all
clc
 
%% Image analysis sample code
 
v = VideoReader('20211122_120358.mp4');
 
v.CurrentTime=2; 

while v.CurrentTime<(v.Duration-12) 
    im0 = readFrame(v); 

    im=im0(525:708,241:1674,:); 
    im=rgb2gray(im); 
    im=imgaussfilt(im,2); 
    im=im2bw(im); 
    im = bwareaopen(im,50); 

    for j=1:size(im,2) 
        h(t,j)=find(im(:,j),1,'last'); 
    end
end
 
ppm = 755; 
Fs = v.FrameRate; 
h_m = h./ppm;
 
% Prepare array
z = h_m(1:2*fix(end/2),1:2*fix(end/2)); 
[m,n] = size(z); 
 
%% FFT sample code
 
% FFT
h_hat=fftshift(fft2(z)); % take fft and shift so that zero frequency is centered
h_hat2=abs(h_hat).^2;
 
% Get appropriate frequencies
f = Fs*((-m/2):(m/2-1))/m; % fft time frequencies
omega = 2*pi*f; % fft time angular frequencies
k = 2*pi*ppm*((-n/2):(n/2-1))/n; % fft wavenumber
 
%% Plots
% Plot the spectrum
[K,OMEGA] = meshgrid(k,omega); % make a meshgrid of the coordinates for plotting
figure(1)
hold on
H=pcolor(K,OMEGA,log10(h_hat2)); % plot a heatmap of log10 of spectrum
set(H,'EdgeColor','none');
colormap('turbo')
colorbar;
caxis([4,15]) % clip colormap limits for easier viewing
xlim([0 300]) % only show 1st quadrant, crop to interesting area
ylim([0 90]) % ylim([0 inf])
xlabel('k')
ylabel('\omega')
hold on
g = 9.81
depth = 0.125
omega_predicted = sqrt(g.*k.*tanh(k.*depth))
plot(k, omega_predicted, 'LineWidth',10,'color','black')
xlim([0 300]) % only show 1st quadrant, crop to interesting area
ylim([0 90])
hold on 
plot(k, sqrt(g.*k),'LineWidth',4,'color','r')
xlim([0 300]) % only show 1st quadrant, crop to interesting area
ylim([0 90])
hold on
plot(k, k.*sqrt(g.*depth),'LineWidth',4,'color','blue')
xlim([0 300]) % only show 1st quadrant, crop to interesting area
ylim([0 90])
title('Wave spectrum')