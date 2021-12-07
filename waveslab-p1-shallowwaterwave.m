close all
clear all
clc
 
%% Image analysis sample code # DEEP WATER WAVE
 
v = VideoReader('20211122_120236.mp4'); 
 
v.CurrentTime=1; 
t=1; 
  
while v.CurrentTime<(v.Duration) 
    im0 = readFrame(v); 
    
    % Process image
    im=im0(525:705,245:1453,:); 
    im=rgb2gray(im); 
    im=imgaussfilt(im,2); 
    im=im2bw(im); 
    im = bwareaopen(im,50); 
    
    for j=1:size(im,2) 
        h(t,j)=find(im(:,j),1,'last'); 
    end
    t=t+1; 
    
end
 
ppm = 755; 
Fs = v.FrameRate; 
h_m = h./ppm;
 
% Prepare array
z = h_m(1:2*fix(end/2),1:2*fix(end/2)); 
[m,n] = size(z); 
 
time_fs = [0:1:m-1]; time = time_fs/Fs;
x_spatial =[0:1:n-1]; x = x_spatial/ppm;
deviation = z-mean(mean(z));
 
%% Propagation speed estimation:
g = 9.81;
depth = mean(mean(z));
c_SW =  sqrt(g.*depth);
 
 
%% Plot
figure(1);
diagram=pcolor(x,time,deviation)
set(diagram,'EdgeColor','none');
colorbar;
 
xlabel('spatial coordinate [m]')
ylabel('time [sec] ')
title('H diagram')
hold on
plot(x,20.2 + c_SW.*x,'LineWidth',3, 'color','r','LineStyle','--')
plot(x,21.2 + c_SW.*x,'LineWidth',3, 'color','r','LineStyle','--')
hold off
