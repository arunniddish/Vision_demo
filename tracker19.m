%%% This code is just a combined version of tracker03.m and tracker05.m
%%% This code will perform real time tracking. The video is live relayed
%%% from webcam. 
%%% The issue i am facing with this code is - all the marker points should
%%% be in the camera view from t = 0 or the first frame.
%%% Live tracking of Caitlin Inchworm for phase delay
%%% This tracking is done on only one point.
%%% Reading and writing data to serial port

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
m = 0;
T_d = 1000;
HIGH = 1;
LOW = 0;
numberOfFrames = 3000;
k = 0;
endpt1 = zeros(numberOfFrames,2);
endpt2 = zeros(numberOfFrames,2);
marker = zeros(numberOfFrames,2);

% Initializing serial port
device = serialport('COM10',9600);

% Initializing camera
cam = webcam;

thisFrame = cam.snapshot;
thisFrame = imresize(thisFrame,[270,480]); % Changing the resolution
newim1 = createMaskpoor(thisFrame);
newim2 = createMaskInchworm_blue(thisFrame);

% Filter out small blobs
newim1 = bwareaopen(newim1, 20);
newim2 = bwareaopen(newim2, 20);

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
  
    % reads the gait cycle number as a char value
    m_char = readline(device);
    
    % converts to double
    m = str2double(m_char);
  
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

%Phase Delay
T_d = round(T_d + K*(marker(m_count-1,1) - marker(m_count-2,1)));
% converts double to string 
Td = num2str(Td);
% writes string to serial port
writeline(device,Td);
end   
     
     
 end
  
    
    
    
    

