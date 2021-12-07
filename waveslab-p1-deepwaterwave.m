close all
clear all

%% Image analysis sample code for DEEP WATER WAVE VIDEO

v = VideoReader('20211122_115959.mp4'); % Make a Matlab VideoReader object from the video file

v.CurrentTime=1; % update current time to skip first 20s of video
t=1; % integer time index for 

while v.CurrentTime<(v.Duration) %loop through frames starting at 10s until 4 seconds before end of video
    % Get image at CurrentTime
    % Each subsequent call of readFrame(v) reads the next frame and updates v.CurrentTime
    im0 = readFrame(v); % original image
    % Process image
    im=im0(590:690,250:1600,:); % crop to an appropriate size
    im=rgb2gray(im); % convert to gray scale
    im=imgaussfilt(im,2); % smoothing to remove noise
    im=im2bw(im); % convert to black and white
    im = bwareaopen(im,50); % remove small glares/bubles by keeping only large objects
    
    % Find free surface
    for j=1:size(im,2) % loop through x positions
        h(t,j)=find(im(:,j),1,'last'); % find the height of the wave as the last element with value 1
    end
    t=t+1; % increment time for entries in h array
end

%% Could reinterpolate h along x-axis based on camera distortion and known 20cm markings in original video frame
ppm = 755
h_m = h./ppm

%% FFT sample code

% Prepare array
z = h_m(1:2*fix(end/2),1:2*fix(end/2)); % make dimensions of array even numbers
[m,n] = size(z); % m is the number of points in the time dimension, n is the number of points in the spatial dimension

% FFT
h_hat=fftshift(fft2(z)); % take fft and shift so that zero frequency is centered
h_hat2=abs(h_hat).^2;

% Get appropriate frequencies
Fs = v.FrameRate; % time sampling rate (frames per second)
f = Fs*((-m/2):(m/2-1))/m; % fft time frequencies
omega = 2*pi*f; % fft time angular frequencies

ppm = 755; % spatial sampling rate (horizontal pixels per meter) (find from known 20cm markings in original video frame)
k = 2*pi*ppm*((-n/2):(n/2-1))/n; % fft wavenumber
%% Plot figure
figure 

[a,b]=size(h_m)

time =[1:1:a]; time_in_s = time/Fs
x_spatial =[0:1:b-1]; x_spatial_to_m = x_spatial/ppm
prop_wave = 28.3 + (1/0.61)*(linspace(0.,1.8,100))
prop_wave_2 = 29.3 + (1/0.61)*(linspace(0.,1.8,100))
hold on 
h_new = pcolor(x_spatial_to_m,time_in_s,(h_m))
set(h_new, 'EdgeColor', 'none')
colorbar;
ylim([800/Fs 1000/Fs]) 
% ylim([0,70])
% xlim([0,1.8])
xlabel('spatial coordinate')
ylabel('time')
hold on
plot(linspace(0.,1.8,100),prop_wave, 'LineWidth',3, 'color','r','LineStyle','--')
plot(linspace(0.,1.8,100),prop_wave_2, 'LineWidth',3, 'color','r','LineStyle','--')

xlabel('spatial coordinate')
ylabel('time')
shading 'flat';