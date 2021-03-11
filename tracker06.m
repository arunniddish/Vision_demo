%%% This code is just a combined version of tracker03.m and tracker05.m
%%% This code will perform real time tracking. The video is live relayed
%%% from webcam. 

%%% The issue i am facing with this code is - all the marker points should
%%% be in the camera view from t = 0 or the first frame.


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
k = 1;
cam = webcam;

folder = pwd; % Current folder path

count = 0;

a = load('cent1.txt');
b = load('cent2.txt');
c = load('cent3.txt');

% Read one frame at a time and find specified color
% for k = 1 : numberOfFrames
while true
    
    if k == 1
		% Enlarge figure to full screen.
		set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
		% Give a name to the title bar.
		set(gcf, 'Name', 'Particle Tracking', 'NumberTitle', 'Off') 
		hCheckbox = uicontrol('Style','checkbox',... 
			'Units', 'Normalized',...
			'String', 'Finish Now',... 
			'Value',0,'Position', [.2 .96 .4 .05], ...
			'FontSize', 14);
        centroid1 = zeros(1000,2);
        centroid2 = zeros(1000,2);
        centroid3 = zeros(1000,2);
	end
    
    thisFrame=cam.snapshot;
    thisFrame = imresize(thisFrame,[270,480]);
    newim = createMaskLR4_1(thisFrame);
    % Filter out small blobs
    newim = bwareaopen(newim, 50);
    % Fill in holes
    newim = imfill(newim, 'holes');
    hImage=subplot(4, 1, 1);
    % Display it.
	imshow(thisFrame);
	axis on;
	caption = sprintf('Original RGB image, frame #%d 0f %d', k, k);
	title(caption, 'FontSize', fontSize);
	drawnow;
    subplot(4,1,2);
    imshow(newim);
    title('Colored Blob Mask', 'FontSize', fontSize);
	drawnow;
    
    [labeledImage, numberOfRegions] = bwlabel(newim);
	if numberOfRegions == 3     % I have changed here to number of regions equal to 3 to avoid errors.
		stats = regionprops(labeledImage, 'BoundingBox', 'Centroid');
        count = count+1;
		% Delete old texts and rectangles
		if exist('hRect', 'var')
			delete(hRect);
		end
		if exist('hText', 'var')
			delete(hText);
		end
		
		% Display the original image again.
		subplot(4, 1, 3); % Switch to original image.
		hImage=subplot(4, 1, 3);
		imshow(thisFrame);
		axis on;
		hold on;
		caption = sprintf('%d blobs found in frame #%d 0f %d', numberOfRegions, k, k);
		title(caption, 'FontSize', fontSize);
		drawnow;
        
        %This is a loop to bound the colored objects in a rectangular box.
		for r = 1 : numberOfRegions
			% Find location for this blob.
			thisBB = stats(r).BoundingBox;
			thisCentroid = stats(r).Centroid;
       
            % Nearest neighbour
            if(count==1 && r==1)
                centroid1(count,:) = thisCentroid;
                filenametowrite = 'cent1.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
            if(count==1 && r==2)
                centroid2(count,:) = thisCentroid;
                 filenametowrite = 'cent2.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
            if(count==1 && r==3)
                centroid3(count,:) = thisCentroid;
                filenametowrite = 'cent3.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
        
            fclose('all');
            
            if(count~=1)
                X1 = [thisCentroid;centroid1(count-1,:)];
                d1 = pdist(X1,'euclidean');
                
                X2 = [thisCentroid;centroid2(count-1,:)];
                d2 = pdist(X2,'euclidean');
            
                X3 = [thisCentroid;centroid3(count-1,:)];
                d3 = pdist(X3,'euclidean');
            
           
        if(d1<d2 && d1<d3)
            centroid1(count,:) = thisCentroid;
                filenametowrite = 'cent1.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d2<d1 && d2<d3)
            centroid2(count,:) = thisCentroid;
                filenametowrite = 'cent2.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
       
        
        if(d3<d1 && d3<d2)
            centroid3(count,:) = thisCentroid;
                filenametowrite = 'cent3.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
            end
       
       
       fclose('all');
       
			hRect(r) = rectangle('Position', thisBB, 'EdgeColor', 'r', 'LineWidth', 2);
			hSpot = plot(thisCentroid(1), thisCentroid(2), 'y+', 'MarkerSize', 10, 'LineWidth', 2)
			hText(r) = text(thisBB(1), thisBB(2)-20, strcat('X: ', num2str(round(thisCentroid(1))), '    Y: ', num2str(round(thisCentroid(2)))));
			set(hText(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
		end
		hold off
		drawnow;
            
        
        
        % Calculating and plotting the body rotation(change in orientation)
        
    if(count~=1) 

    % Reading previous coordinate of the points
       vi1x = a(1,1);
       vi1y = a(1,2);
       vi2x = b(1,1);
       vi2y = b(1,2);
       vi3x = c(1,1);
       vi3y = c(1,2);

    % Reading the current coordinate of the points
      vf1x = a(count,1);
      vf1y = a(count,2);
      vf2x = b(count,1);
      vf2y = b(count,2);
      vf3x = c(count,1);
      vf3y = c(count,2);
    

     % Angle between 2 vectors - initial state and final state
        theta_1 = orientation(vi2x,vi1x,vi2y,vi1y,vf2x,vf1x,vf2y,vf1y);
        
        theta_2 = orientation(vi3x,vi1x,vi3y,vi1y,vf3x,vf1x,vf3y,vf1y);

        theta_3 = orientation(vi3x,vi2x,vi3y,vi2y,vf3x,vf2x,vf3y,vf2y);



    %  Mean of all theta's
        theta = ((theta_1 + theta_2 + theta_3)/3);
    
        %Plot    
        subplot(4,1,4);
        plot(count,theta,'r*');
            xlim([0 260]);
            ylim([-80 200]);
            drawnow 
            hold on 
    end
    end
    k = k+1;
   
	% See if they want to bail out
	if get(hCheckbox, 'Value')
		% Finish now checkbox is checked.
		msgbox('Done with processing.');
		return;
	end
end
msgbox('Done with processing.');

%%
% Plots of the tracked points
a = load('cent1.txt');
b = load('cent2.txt');
c = load('cent3.txt');
figure(2)
plot(a(:,1));
hold on
plot(a(:,2));
legend('Point1 x','Point1 y');
figure(3)
plot(b(:,1));
hold on
plot(b(:,2));
legend('Point2 x','Point2 y');
figure(4)
plot(c(:,1));
hold on
plot(c(:,2));
legend('Point3 x','Point3 y');

%%

% Reading initial coordinate of the points
vi1x = a(1,1);
vi1y = a(1,2);
vi2x = b(1,1);
vi2y = b(1,2);
vi3x = c(1,1);
vi3y = c(1,2);

% Reading the final coordinate of the points
vf1x = a(k,1);
vf1y = a(k,2);
vf2x = b(k,1);
vf2y = b(k,2);
vf3x = c(k,1);
vf3y = c(k,2);

% Vector v12 ( Arbitary vector pointing from 1 to 2)
vi12 = [(vi2x-vi1x) (vi2y-vi1y) 0];
vf12 = [(vf2x-vf1x) (vf2y-vf1y) 0];

% Angle between 2 vectors
theta_1 = atan2(norm(cross(vi12,vf12)),dot(vi12,vf12));


%Vector v13 ( Arbitary vector pointing from 1 to 3)
vi13 = [(vi3x-vi1x) (vi3y-vi1y) 0];  
vf13 = [(vf3x-vf1x) (vf3y-vf1y) 0];
% I am commenting this part because the deviation from other two angle is
% higher. This increases the error percentage
% % Angle between 2 vectors
% theta_2 = atan2(norm(cross(vi13,vf12)),dot(vi13,vf12));

% Vector v23 ( Arbitary vector pointing from 2 to 3)
vi23 = [(vi3x-vi2x) (vi3y-vi2y) 0];
vf23 = [(vf3x-vf2x) (vf3y-vf2y) 0];

% Angle between 2 vectors
theta_3 = atan2(norm(cross(vi23,vf23)),dot(vi23,vf23));

% Mean of all theta's
theta = ((theta_1 + theta_3)/2)*180/pi;

% Identifying the respective points
gamma_1 = (atan2(norm(cross(vf12,vf13)),dot(vf12,vf13)))*180/pi;
gamma_2 = (atan2(norm(cross(vf12,vf23)),dot(vf12,vf23)))*180/pi;
% gamma_3 = (atan2(norm(cross(vf13,vf23)),dot(vi13,vf23)))*180/pi;
gamma_3 = 180-(gamma_1 + gamma_2);


if((20<gamma_1)&&(gamma_1<35))
filenametowrite = 'cent1.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
fprintf(fileid,'Point 1');
end
 
if((80<gamma_2)&&(gamma_2<100)) 
filenametowrite = 'cent2.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
fprintf(fileid,'Point 2');
end
 
if((50<gamma_3)&&(gamma_3<70))
filenametowrite = 'cent3.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
fprintf(fileid,'Point 3');
end

% Calculating the Rotation matrix
Ai = [vi12;vi13;vi23];
Af = [vf12;vf13;vf23];
Rot_m = pinv(Af)*Ai;











