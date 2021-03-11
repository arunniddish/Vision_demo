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

% Specify input video file name.
folder = pwd; % Current folder path
baseFileName = 'track180.mp4';  % Video to track
fullFileName = fullfile(folder, baseFileName);

% Check if the video file actually exists in the current folder or on the search path.
if ~exist(fullFileName, 'file')
	% File doesn't exist -- didn't find it there.  Check the search path for it.
	fullFileNameOnSearchPath = baseFileName; % No path this time.
	if ~exist(fullFileNameOnSearchPath, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end

% Instantiate a video reader object for this video.
videoObject = VideoReader(fullFileName);

% Setup other parameters
numberOfFrames = videoObject.NumberOfFrame;

% Capture one of the frame to evaluate color profile

% b = read(videoObject,8);
% imwrite(b,'capture1.jpg');
% 
% pause();

% Read one frame at a time and find specified color
for k = 1 : numberOfFrames
    
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
        centroid1 = zeros(numberOfFrames,2);
        centroid2 = zeros(numberOfFrames,2);
        centroid3 = zeros(numberOfFrames,2);
	end
    
    
    thisFrame=read(videoObject,k);
    newim = createMask(thisFrame);
    % Filter out small blobs
    newim = bwareaopen(newim, 200);
    % Fill in holes
    newim = imfill(newim, 'holes');
    hImage=subplot(4, 1, 1);
    % Display it.
	imshow(thisFrame);
	axis on;
	caption = sprintf('Original RGB image, frame #%d 0f %d', k, numberOfFrames);
	title(caption, 'FontSize', fontSize);
	drawnow;
    subplot(4,1,2);
    imshow(newim);
    title('Colored Blob Mask', 'FontSize', fontSize);
	drawnow;
    
    [labeledImage, numberOfRegions] = bwlabel(newim);
	if numberOfRegions >= 1
		stats = regionprops(labeledImage, 'BoundingBox', 'Centroid');
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
		caption = sprintf('%d blobs found in frame #%d 0f %d', numberOfRegions, k, numberOfFrames);
		title(caption, 'FontSize', fontSize);
		drawnow;
        
        %This is a loop to bound the colored objects in a rectangular box.
		for r = 1 : numberOfRegions
			% Find location for this blob.
			thisBB = stats(r).BoundingBox;
			thisCentroid = stats(r).Centroid;
            
            % Nearest neighbour
            if(k==1 && r==1)
                centroid1(k,:) = thisCentroid;
                filenametowrite = 'cent1.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
            if(k==1 && r==2)
                centroid2(k,:) = thisCentroid;
                 filenametowrite = 'cent2.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
            if(k==1 && r==3)
                centroid3(k,:) = thisCentroid;
                filenametowrite = 'cent3.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
        
            if(k~=1)
                X1 = [thisCentroid;centroid1(k-1,:)];
                d1 = pdist(X1,'euclidean');
                
                X2 = [thisCentroid;centroid2(k-1,:)];
                d2 = pdist(X2,'euclidean');
            
                X3 = [thisCentroid;centroid3(k-1,:)];
                d3 = pdist(X3,'euclidean');
            
           
        if(d1<d2 && d1<d3)
            centroid1(k,:) = thisCentroid;
                filenametowrite = 'cent1.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d2<d1 && d2<d3)
            centroid2(k,:) = thisCentroid;
                filenametowrite = 'cent2.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
       
        
        if(d3<d1 && d3<d2)
            centroid3(k,:) = thisCentroid;
                filenametowrite = 'cent3.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
            end
            
			hRect(r) = rectangle('Position', thisBB, 'EdgeColor', 'r', 'LineWidth', 2);
			hSpot = plot(thisCentroid(1), thisCentroid(2), 'y+', 'MarkerSize', 8, 'LineWidth', 2);
			hText(r) = text(thisBB(1), thisBB(2)-20, strcat('X: ', num2str(round(thisCentroid(1))), '    Y: ', num2str(round(thisCentroid(2)))));
			set(hText(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 8, 'Color', 'yellow');
		end
		hold off
		drawnow;
       
        
    % Calculating and plotting the body rotation(change in orientation)
        
    if(k~=1 && r==3)  
        a = load('cent1.txt');
        b = load('cent2.txt');
        c = load('cent3.txt');
    % Reading previous coordinate of the points
       vi1x = a(1,1);
       vi1y = a(1,2);
       vi2x = b(1,1);
       vi2y = b(1,2);
       vi3x = c(1,1);
       vi3y = c(1,2);

    % Reading the current coordinate of the points
      vf1x = a(k,1);
      vf1y = a(k,2);
      vf2x = b(k,1);
      vf2y = b(k,2);
      vf3x = c(k,1);
      vf3y = c(k,2);
    

     % Angle between 2 vectors - initial state and final state
        theta_1 = orientation(vi2x,vi1x,vi2y,vi1y,vf2x,vf1x,vf2y,vf1y);
        
        theta_2 = orientation(vi3x,vi1x,vi3y,vi1y,vf3x,vf1x,vf3y,vf1y);

        theta_3 = orientation(vi3x,vi2x,vi3y,vi2y,vf3x,vf2x,vf3y,vf2y);



    %  Mean of all theta's
        theta = ((theta_1 + theta_2 + theta_3)/3);
    
        %Plot    
        subplot(4,1,4);
        plot(k,theta,'r*');
            xlim([0 260]);
            ylim([-180 200]);
            drawnow 
            hold on 
    end
    end
    
    
	% See if you want to bail out
	if get(hCheckbox, 'Value')
		% Finish now checkbox is checked.
		msgbox('Done with processing.');
		return;
	end
end
msgbox('Done with processing.');


% Plots of the tracked points
% a = load('cent1.txt');
% b = load('cent2.txt');
% c = load('cent3.txt');
figure(2)
subplot(3,1,1);
plot(a(:,1));
hold on
plot(a(:,2));
legend('Point1 x','Point1 y');
subplot(3,1,2);
plot(b(:,1));
hold on
plot(b(:,2));
legend('Point2 x','Point2 y');
subplot(3,1,3);
plot(c(:,1));
hold on
plot(c(:,2));
legend('Point3 x','Point3 y');



% Reading initial coordinate of the points
% vi1x = a(1,1);
% vi1y = a(1,2);
% vi2x = b(1,1);
% vi2y = b(1,2);
% vi3x = c(1,1);
% vi3y = c(1,2);

% Reading the final coordinate of the points
% vf1x = a(k,1);
% vf1y = a(k,2);
% vf2x = b(k,1);
% vf2y = b(k,2);
% vf3x = c(k,1);
% vf3y = c(k,2);



% Vector v12 ( Arbitary vector pointing from 1 to 2)
% vi12 = [(vi2x-vi1x) (vi2y-vi1y) 0];
% vi12 = vi12/norm(vi12);
vf12 = [(vf2x-vf1x) (vf2y-vf1y) 0];
% vf12 = vf12/norm(vf12);
% % [VV] - Normalize vi12, vf12
% % Angle between 2 vectors - initial state and final state
% theta_1 = atan2(norm(cross(vi12,vf12)),dot(vi12,vf12));
% theta_1 = radtodeg(theta_1);



%Vector v13 ( Arbitary vector pointing from 1 to 3)
% vi13 = [(vi3x-vi1x) (vi3y-vi1y) 0];
% vi13 = vi13/norm(vi13);
vf13 = [(vf3x-vf1x) (vf3y-vf1y) 0];
% vf13 = vf13/norm(vf13);
% % [VV] - Normalize vi12, vf12
% % Angle between 2 vectors - initial state and final state
% theta_2 = atan2(norm(cross(vi13,vf13)),dot(vi13,vf13));
% theta_2 = radtodeg(theta_2);



% Vector v23 ( Arbitary vector pointing from 2 to 3)
% vi23 = [(vi3x-vi2x) (vi3y-vi2y) 0];
% vi23 = vi23/norm(vi23);
vf23 = [(vf3x-vf2x) (vf3y-vf2y) 0];
% vf23 = vf23/norm(vf23);
% % Angle between 2 vectors - initial state and final state
% theta_3 = atan2(norm(cross(vi23,vf23)),dot(vi23,vf23));
% theta_3 = radtodeg(theta_3);



% Mean of all theta's
% theta = ((theta_1 + theta_2 + theta_3)/3);



% Identifying the respective points
gamma_1 = radtodeg((atan2(norm(cross(vf12,vf13)),dot(vf12,vf13))));
gamma_2 = radtodeg((atan2(norm(cross(vf12,vf23)),dot(vf12,vf23))));
gamma_3 = radtodeg((atan2(norm(cross(vf13,vf23)),dot(vf13,vf23))));



% Specifying the conditions of the sides of the triangle to retain the
% respective points
figure(3)
imshow(thisFrame);
title('Retaining marker points', 'FontSize', fontSize);
hold on

if((20<gamma_1)&&(gamma_1<30))
filenametowrite = 'cent1.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
fprintf(fileid,'Point A');
plot(a(k,1),a(k,2),'y+','MarkerSize',4,'LineWidth',1);
text(a(k,1),a(k,2),'A','Color','black','FontSize',10);
end


if((85<gamma_2)&&(gamma_2<95)) 
filenametowrite = 'cent2.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
fprintf(fileid,'Point B');
plot(b(k,1),b(k,2),'y+','MarkerSize',4,'LineWidth',1);
text(b(k,1),b(k,2),'B','Color','black','FontSize',10);
end


if((60<gamma_3)&&(gamma_3<70))
filenametowrite = 'cent3.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
fprintf(fileid,'Point C');
plot(c(k,1),c(k,2),'y+','MarkerSize',4,'LineWidth',1);
text(c(k,1),c(k,2),'C','Color','black','FontSize',10);
end

% Calculating the Rotation matrix
Ai = [vi12;vi13;vi23];
Af = [vf12;vf13;vf23];
Rot_m = pinv(Af)*Ai;


% Plotting the change in orientation













