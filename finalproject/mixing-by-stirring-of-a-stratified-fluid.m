% clc
% clf
% close all
% clear all

%% Context

% ------------------------------------------------------------

% Project mates: Aurora Leeson, Scott Martin, and Arief Suryo
% Code written for OCEAN511 Final Project 

% --------------- Project description ------------------------

% Mixing by stirring of a stratified fluid (GFD lab). Explore 
% the dependence of the turbulent tracer diffusivity on the 
% stratification of the fluid. Quantify the mixing efficiency 
% for different stratifications.  

% Not the final version of the code. Some mistakes here.
% ------------------------------------------------------------

%% Image analysis sample code (borrowed code from wave_sample_code.m)

v = VideoReader('IMG_6493_MOV_SparkVideo.mp4'); % Make a Matlab VideoReader object from the video file

v.CurrentTime=16; % update current time to skip first 16s of video
t=1; % integer time index for 

fig = figure(); 

while v.CurrentTime<(v.Duration-4) %loop through frames starting at 16s until 4 seconds before end of video
    % Get image at CurrentTime
    % Each subsequent call of readFrame(v) reads the next frame and updates v.CurrentTime
    im0 = readFrame(v); % original image
    
    % Process image
     im=rgb2gray(im0); % convert to gray scale
     im = imrotate(im, -2.5); % straigten image
     
%      keyboard()
     
    if 1 == 1 %mod(t,200) == 1                                              
        
        % Crop to each container
%         im1 = im(340:620,100:480,:);
%         im2 = im(340:620,470:850,:);
%         im3 = im(340:620,1170:1550,:);
%         im4 = im(340:620,1540:1920,:);

        % Crop to small region in each contain to get pixel value
        im1 = im(340:620,200:250,:);
        im2 = im(340:620,570:620,:);
        im3 = im(340:620,1270:1320,:);
        im4 = im(340:620,1640:1690,:);

        %Option for plotting
%         subplot(1,4,1)
%         image(im1)
%         colormap('gray')
%         title({['t =' num2str(round(t,2))], 'rho = 1.02 g/mL'})
%     
%         subplot(1,4,2)
%         image(im2)
%         colormap('gray')
%         title({['T=' num2str(round(v.CurrentTime,2)) 's'], 'rho = 1.06 g/mL'})
%     
%         subplot(1,4,3)
%         image(im3)
%         colormap('gray')
%         title({['T=' num2str(round(v.CurrentTime,2)) 's'], 'rho = 1.13 g/mL'})
%     
%         subplot(1,4,4)
%         image(im4)
%         colormap('gray')
%         title({['T=' num2str(round(v.CurrentTime,2)) 's'], 'rho = 1.16 g/mL'})
%     
%         pause(1/v.FrameRate); % pause before showing next frame
%         set(gcf,'Color','White')

        % Get initial data
        if t == 1
            % get average value in to build gradient in y
            px1_t0 = mean(im1,2);
            px2_t0 = mean(im2,2);
            px3_t0 = mean(im3,2);
            px4_t0 = mean(im4,2);
    
        % Get data for some arbitrary time
        elseif t == 801
            % get average value in pixel intensity to build gradient
            px1_t1 = mean(im1,2);
            px2_t1 = mean(im2,2);
            px3_t1 = mean(im3,2);
            px4_t1 = mean(im4,2);

        % Get data for some arbitrary time
        elseif t == 2201
            % get average value in pixel intensity to build gradient
            px1_t2 = mean(im1,2);
            px2_t2 = mean(im2,2);
            px3_t2 = mean(im3,2);
            px4_t2 = mean(im4,2);

        % Get data for some arbitrary time
        elseif t == 4001
            % get average value in pixel intensity to build gradient
            px1_t3 = mean(im1,2);
            px2_t3 = mean(im2,2);
            px3_t3 = mean(im3,2);
            px4_t3 = mean(im4,2);

        % Get data for some arbitrary time
        elseif t == 8001
            % get average value in pixel intensity to build gradient
            px1_t4 = mean(im1,2);
            px2_t4 = mean(im2,2);
            px3_t4 = mean(im3,2);
            px4_t4 = mean(im4,2);

        % Get data for some arbitrary time
        elseif t == 10001
            % get average value in pixel intensity to build gradient
            px1_t5 = mean(im1,2);
            px2_t5 = mean(im2,2);
            px3_t5 = mean(im3,2);
            px4_t5 = mean(im4,2);

        end

        px1(t) = mean(im1(25,:));
        px2(t) = mean(im2(25,:));
        px3(t) = mean(im3(25,:));
        px4(t) = mean(im4(25,:));

        time(t) = t;

        clf(fig); % clear figure
    end
    
    t=t+1; % increment time for entries in h array
    if mod(t,200) == 1
        disp(t)
    end
    
end

depth = 1:1:length(px1_t0); % depth in pixels

%% Making some figures
figure(2)
clf

subplot(1,4,1)
plot(px1_t0,-1*depth)
hold on
plot(px1_t1,-1*depth)
plot(px1_t2,-1*depth)
plot(px1_t3,-1*depth)
plot(px1_t4,-1*depth)
plot(px1_t5,-1*depth)
legend('t = 0', ['t = ' num2str(round(801/v.FrameRate))], ...
    ['t = ' num2str(round(2201/v.FrameRate))], ...
    ['t = ' num2str(round(4001/v.FrameRate))], ...
    ['t = ' num2str(round(8001/v.FrameRate))], ...
    ['t = ' num2str(round(10001/v.FrameRate))], 'Location', 'southeast')
xlabel('Pixel Intensity')
ylabel('Depth (pixels)')
title('rho = 1.02 g/mL')

subplot(1,4,2)
plot(px2_t0,-1*depth)
hold on
plot(px2_t1,-1*depth)
plot(px2_t2,-1*depth)
plot(px2_t3,-1*depth)
plot(px2_t4,-1*depth)
plot(px2_t5,-1*depth)
legend('t = 0', ['t = ' num2str(round(801/v.FrameRate))], ...
    ['t = ' num2str(round(2201/v.FrameRate))], ...
    ['t = ' num2str(round(4001/v.FrameRate))], ...
    ['t = ' num2str(round(8001/v.FrameRate))], ...
    ['t = ' num2str(round(10001/v.FrameRate))], 'Location', 'southeast')
xlabel('Pixel Intensity')
ylabel('Depth (pixels)')
title('rho = 1.06 g/mL')

subplot(1,4,3)
plot(px3_t0,-1*depth)
hold on
plot(px3_t1,-1*depth)
plot(px3_t2,-1*depth)
plot(px3_t3,-1*depth)
plot(px3_t4,-1*depth)
plot(px3_t5,-1*depth)
legend('t = 0', ['t = ' num2str(round(801/v.FrameRate))], ...
    ['t = ' num2str(round(2201/v.FrameRate))], ...
    ['t = ' num2str(round(4001/v.FrameRate))], ...
    ['t = ' num2str(round(8001/v.FrameRate))], ...
    ['t = ' num2str(round(10001/v.FrameRate))], 'Location', 'southeast')
xlabel('Pixel Intensity')
ylabel('Depth (pixels)')
title('rho = 1.13 g/mL')

subplot(1,4,4)
plot(px4_t0,-1*depth)
hold on
plot(px4_t1,-1*depth)
plot(px4_t2,-1*depth)
plot(px4_t3,-1*depth)
plot(px4_t4,-1*depth)
plot(px4_t5,-1*depth)
legend('t = 0', ['t = ' num2str(round(801/v.FrameRate))], ...
    ['t = ' num2str(round(2201/v.FrameRate))], ...
    ['t = ' num2str(round(4001/v.FrameRate))], ...
    ['t = ' num2str(round(8001/v.FrameRate))], ...
    ['t = ' num2str(round(10001/v.FrameRate))], 'Location', 'southeast')
xlabel('Pixel Intensity')
ylabel('Depth (pixels)')
title('rho = 1.16 g/mL')

set(gcf,'Color','White')


beaker1 = [px1_t0(25), px1_t1(25), px1_t2(25), px1_t3(25), px1_t4(25), px1_t5(25)];
beaker2 = [px2_t0(25), px2_t1(25), px2_t2(25), px2_t3(25), px2_t4(25), px2_t5(25)];
beaker3 = [px3_t0(25), px3_t1(25), px3_t2(25), px3_t3(25), px3_t4(25), px3_t5(25)];
beaker4 = [px4_t0(25), px4_t1(25), px4_t2(25), px4_t3(25), px4_t4(25), px4_t5(25)];
%time = [0, 801/v.FrameRate, 2201/v.FrameRate 4001/v.FrameRate 8001/v.FrameRate 10001/v.FrameRate];

figure(3);
clf
% plot(time, beaker1)
% hold on
% plot(time,beaker2)
% plot(time,beaker3)
% plot(time,beaker4)
px1 = px1 - px1(end);
plot(time/v.FrameRate, px1/px1(1))
hold on
px2 = px2 - px2(end);
plot(time/v.FrameRate, px2/px2(1))
px3 = px3 - px3(end);
plot(time/v.FrameRate, px3/px3(1))
px4 = px4 - px4(end);
plot(time/v.FrameRate, px4/px4(1))
% ylim([0 500]);
legend('Beaker 1','Beaker 2', 'Beaker 3', 'Beaker 4')
xlabel('time(s)')
ylabel('normalized pixel intensity')
set(gcf,'Color','White')

%% Modified diffusion code (provided from class files)
N=20; % number of grid boxes

L=9.5e-2; % length of the domain
dx=(L/N); x=dx/2:dx:L-dx/2;

T0 = ones([N,1]);
T0((N/2+1):N) = -1;

T=T0;

D=0.000021; % diffusivity coefficient

Nt=10000; 

% Courant stability requires 0.5 prefactor minimum #no need for this
dt=0.25*dx^2/D;
disp('time step:')
disp(dt)
time=0;

T_top=zeros(1,Nt);
time_series = zeros(1,Nt);

Txx=zeros(size(T));

fig=figure; 
i=1;
while T(1)-T(end)>0.01        %for i=0:Nt
%for i =1:Nt    
    %plotting 
    if mod(i,10)==0
        clf(fig)
        plot(x,T0,'k-'); hold on
        plot(x,T);
        %plot(time,T(1),'x'); hold on;
        legend('Initial','Simulation')
         pause(0.1)
    end
    
    % calculating Txx using a simple central difference scheme, O(dx^2)
    Txx(2:end-1)=(T(3:end)-2*T(2:end-1) + T(1:end-2))/dx^2;
    %Txx(1)=(T(2)-2*T(1) + T(end))/dx^2;
    %Txx(end)=(T(1)-2*T(end) + T(end-1))/dx^2;

    %no0flux 
    Txx(1)=(T(2)-2*T(1) + T(1))/dx^2;
    Txx(end)=(T(end)-2*T(end) + T(end-1))/dx^2;

    % updating temperature (time-stepping)
    T=T+D*Txx*dt;    
    T_top(i)=T(1);
    time_series(i)=time;

    time=i*dt;
    i=i+1;  
end
disp('time to mix:')
disp(time)
%% Making some figures
T_top = T_top(1:i);
time_series = time_series(1:i);

figure;
plot(T_top)
%% Useful arrays
% delta GPE = 
% N2 = (-g/mean_dens)*(dens_distilled_water - dens_salty_part)/(2H)

N2 = -(g/mean_rho)*((1000-rho_salt)/H);
gpe = -(1/8)*g*(rho_salt_1000)*H^2;

time_in_sec = time/v.FrameRate
b_val = [-0.0208 -0.0155 -0.0109 -0.0057]
% N_squared = [4.088 12.027 25.027 30.586]
% gpe = [-0.22 -0.66 -1.44 -1.77]
mean_dens = [1010 1030 1065 1080]
densest_dens = [1020 1060 1130 1160]

% delta_t = [120.9400 169.5000 229.500 437.2500]
K = [0.21e-04 0.1570e-04 0.11e-04 0.0583e-04]

%(1/1-(1/e))GPE * b = dGPE / dT


% mixing happens in interface
% nonlinear color scheme --> density --> overestimating Kp
% time that you pick up --> biasing
%% Make some figures
figure;
y_axis = mean_dens.*N_squared.*K
x_axis = gpe./delta_t

plot(N_squared.*K,gpe./delta_t)
%% Delete one and try this
figure;
plot(N_squared(1:3).*K(1:3),gpe(1:3)./delta_t(1:3))
x_axis_first3 = N_squared(1:3).*K(1:3)
y_axis_first3 = gpe(1:3)./delta_t(1:3)
