%%% This code is just a combined version of tracker03.m and tracker05.m
%%% This code will perform real time tracking. The video is live relayed
%%% from webcam. 
%%% The issue i am facing with this code is - all the marker points should
%%% be in the camera view from t = 0 or the first frame.
%%% Live tracking of Caitlin Inchworm for phase delay
%%% This tracking is done on only one point.

% Demo to track any color based on 'createMask' function.  Finds and annotates centroid and bounding box of any colored blobs.
% Modify thresholds to detect different colors. In this case we are finding
% magenta colors and create blobs around it.


clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;
folder = pwd;
count = 0;
m_count = 0;
HIGH = 1;
LOW = 0;
numberOfFrames = 3000;
endpt1 = zeros(numberOfFrames,2);
endpt2 = zeros(numberOfFrames,2);
marker = zeros(numberOfFrames,2);

% Initializing arduino
a = arduino()

% Arduino data initialization

% define Seeduino pin numbers. These can be reversed to change the locomotion direction
 left = 8;
 right = 7;
 left_led = 9;
 right_led = 10;

% boolean flags recording whether LEDs are currently on
left_led_on = false;
right_led_on = false;

% max pulse width modulation value (out of 225). This results in a max duty cycle of 60/225 = 27%
 pwm_val_l =75;   % left actuator pwm
 pwm_val_r = 85;    % right actuator pwm

% actuation time parameters
rise_time = 3000; %  T_r in milliseconds
pulse_width = 1000; % T_w in milliseconds, time at max actuation
time_delay = 1500; % T_d in milliseconds, time between left and right actuation
time_period = 10000; % T in milliseconds, total time for actuation waveform

% number of steps for data output and analog value writing
time_step = 10; % milliseconds
N_r = rise_time/ time_step;
N = time_period / time_step; 
N_d = time_delay / time_step;
N_w = pulse_width / time_step;

% initialize actuators as being off (0 PWM)
left_signal = 0;
right_signal = 0;

i = 0;
% gait (cycle) number
m = 0;  % To modify
% current values in mA
left_current = 0;
right_current = 0;

% declare functions 
% void set_actuation_signal();
% void record_actuation_signal();
% void run_calibration_routine();

pause(2000);




% Initializing camera
cam = webcam;

thisFrame = cam.snapshot;
thisFrame = imresize(thisFrame,[270,480]); % Changing the resolution
newim1 = createMaskInchworm2(thisFrame);
newim2 = createMaskInchworm_blue(thisFrame);

% Filter out small blobs
newim1 = bwareaopen(newim1, 10);
newim2 = bwareaopen(newim2, 10);

% Fill in holes
newim1 = imfill(newim1, 'holes');
newim2 = imfill(newim2, 'holes');

hImage=subplot(3, 1, 1);
% Display it.
imshow(thisFrame);
axis on;
caption = sprintf('Original image, frame #%d 0f %d', k, numberOfFrames);
title(caption, 'FontSize', fontSize);
drawnow;

[labeledImage, numberOfRegions] = bwlabel(newim1);
if numberOfRegions == 2
    count = count+1;
    stats1 = regionprops(labeledImage, 'BoundingBox', 'Centroid');
    
		% Delete old texts and rectangles
		if exist('hRect', 'var')
			delete(hRect);
		end
		if exist('hText', 'var')
			delete(hText);
        end
end
for r1 = 1 : numberOfRegions 
thisBB = stats1(r1).BoundingBox;
thisCentroid = stats1(r1).Centroid;
thisCentroid = [thisCentroid(1) thisCentroid(2)]*pix2mm;

            % Nearest neighbour
            if(r1==1)
                endpt1(count,:) = thisCentroid;
            end
            
            if(r1==2)
                endpt2(count,:) = thisCentroid;
            end
end

                XX = [endpt1(count,1);endpt2(count,1)];
                dx = pdist(XX,'euclidean');
                
                %Pixel to mm conversion
                pix2mm = 89/dx;
    
    
[labeledImage, numberOfRegions] = bwlabel(newim2);
if numberOfRegions == 1
    m_count = m_count + 1;
   stats2 = regionprops(labeledImage, 'BoundingBox', 'Centroid');
    
		% Delete old texts and rectangles
		if exist('hRect', 'var')
			delete(hRect);
		end
		if exist('hText', 'var')
			delete(hText);
        end
end 
    
for r2 = 1 : numberOfRegions 
thisBB = stats2(r2).BoundingBox;
thisCentroid = stats2(r2).Centroid;
thisCentroid = [thisCentroid(1) thisCentroid(2)]*pix2mm;

            % Nearest neighbour
            if(r2==1)
                marker(m_count,:) = thisCentroid;
            end
          
end

    
while(true)
    %% Put Arduino code here
    
    % each loop is one locomotion gait cycle. It can be divided into four actuation blocks: left limb (2s), both limbs (2s), right limb(2s), no limbs (6s)
for i = 1:N % each loop (cycle) runs for N steps (T ms)
% First part of loop: sets the current for left and right limbs
	  % start with left limb actuation
	if(i <= N_r) 
	  left_signal = i*pwm_val_l/N_r; % if time <= rise time, do ramp input for left limb up to pwm_val_l
    end

	if(i > N_d && i <= N_d+N_r) % if time> time_delay, do ramp input for right limb up to pwm_val_r for duration of rise_time 
	  right_signal = (i-N_d)*pwm_val_r/N_r;
    end

	if(i == N_r+1) 
	  left_signal = pwm_val_l; % as soon as rise_time has passed, set left limb to max current pwm_val_l
		% statements
    end
	if (i == N_d+N_r+1) % if time_delay + rise_time has passed (right limb has reached max pwm_val_r), hold at max current pwm_val_r
	  right_signal = pwm_val_r;
		% statements
    end
	if(i == N_r+N_w+1) % if rise_time+pulse_width has passed (left limb has completed actuation), send left limb to 0 
	  left_signal = 0;
		% statements
    end
	if (i == N_d + N_r+N_w+1) % if rise_time+pulse_width +time_delay has passed (right limb has completed actuation), send right limb to 0 
	  right_signal = 0;
    end


%Second part of loop: actual functions to write current to the outputs and record in serial monitor / plotter
%   set_actuation_signal();

  writePWMVoltage(a,left,left_signal);
  writePWMVoltage(a,right,right_signal);


%   record_actuation_signal(); % This function includes the time step 
  
%    delay(time_step);
pause(time_step);
   left_current =  left_signal*242/pwm_val_l;
   right_current = right_signal*247/pwm_val_r;

% Third part of loop:control LEDs to signal max actuation
	  if(left_led_on == false && left_signal == pwm_val_l) 
		writeDigitalPin(a,left_led, HIGH);
		left_led_on = true;
      end
	  if(right_led_on == false && right_signal == pwm_val_r) 
		writeDigitalPin(a,right_led, HIGH);
		right_led_on = true;
      end
	  if(left_led_on == true && left_signal == 0)
		writeDigitalPin(a,left_led, LOW);
		left_led_on = false;
      end
	  if(right_led_on == true && right_signal == 0) 
		writeDigitalPin(a,right_led, LOW);
		right_led_on = false;
      end
end
m = m+1;

% Function Definitions
% void set_actuation_signal() 
  % This function sets the left limb output to left_signal (PWM) and the right limb output to right_signal (PWM) 
  % This is considered an analog output, but is really just a pulse-modulated digital output
%   writePWMVoltage(a,left,left_signal);
%   writePWMVoltage(a,right,right_signal);



% void record_actuation_signal() 
  % This function sends the left_signal and right_signal PWM values to the Serial Monitor and Serial Plotter in increments of 
  % time_step milliseconds 
%    delay(time_step);
%    left_current = (float) left_signal*242/pwm_val_l;
%    right_current = (float) right_signal*247/pwm_val_r;
% SerialUSB.print("Left_Current:");
% SerialUSB.print(left_current);
% SerialUSB.print(" , ");
% SerialUSB.print("Right_Current:");
% SerialUSB.print(right_current);
% SerialUSB.print(" , ");
% SerialUSB.println(300);



    
%% 
    
    if(m>=2)
        
        [labeledImage, numberOfRegions] = bwlabel(newim2);
        if numberOfRegions == 1
           m_count = m_count + 1;
           stats2 = regionprops(labeledImage, 'BoundingBox', 'Centroid');
    
		% Delete old texts and rectangles
		if exist('hRect', 'var')
			delete(hRect);
		end
		if exist('hText', 'var')
			delete(hText);
        end
        end
    

        for r2 = 1 : numberOfRegions 
thisBB = stats1(r2).BoundingBox;
thisCentroid = stats1(r2).Centroid;
thisCentroid = [thisCentroid(1)-x_correction thisCentroid(2)]*pix2mm;

            % Nearest neighbour
            if(r2==1)
                marker(m_count,:) = thisCentroid;
            end
          
end

% Update T_d

T_d = T_d + K*(marker(m_count-1,1) - marker(m_count-2,1));
    end
end
    
    HIGH = 1;
    LOW = 0;
    while(true)
       
        writeDigitalPin(a,'D7', 1);
        writeDigitalPin(a,'D8', 1);
%         pause(10);
%         writeDigitalPin(a,'D7', 1);
%         writeDigitalPin(a,'D8', 1);
%     pause(10);
    end
    
 %%
 a = arduino('COM11','Uno')
 device = serialport('COM10',9600);   
 
 write(device,'2','uint8')  
 
    
    
    
    
    
    
    
    
    


