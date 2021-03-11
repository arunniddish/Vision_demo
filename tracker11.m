%%% This code is same as tracker 07
%%% This code is to track Caitlin's video - Inch Worm
%%% The Inch worm has 7 markers on it. Therefore 7 different markers are
%%% tracked seperately. 
%%% Right only tendon actuation


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
baseFileName = 'Inchworm_3_Right_only.mp4';  % Video to track
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

count = 0;

% Correction factor/values - x axis and y axis

x_correction = 0;
y_correction = 270.5;

% Pixel to mm conversion

pix2mm = 1;  % This is only for this particular video. Should be updated for different videos



% Read one frame at a time and find specified color
for k = 314
    
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
        centroid4 = zeros(numberOfFrames,2);
        centroid5 = zeros(numberOfFrames,2);
        centroid6 = zeros(numberOfFrames,2);
        centroid7 = zeros(numberOfFrames,2);
        
        a = load('cent1.txt');
        b = load('cent2.txt');
        c = load('cent3.txt');
        d = load('cent4.txt');
        e = load('cent5.txt');
        f = load('cent6.txt');
        g = load('cent7.txt');
        
	end
    
    
    thisFrame=read(videoObject,k);
    thisFrame = imresize(thisFrame,[270,480]);
    newim = createMaskInchworm2(thisFrame);
    % Filter out small blobs
    newim = bwareaopen(newim, 10);
    % Fill in holes
    newim = imfill(newim, 'holes');
    hImage=subplot(3, 1, 1);
    % Display it.
	imshow(thisFrame);
	axis on;
	caption = sprintf('Original RGB image, frame #%d 0f %d', k, numberOfFrames);
	title(caption, 'FontSize', fontSize);
	drawnow;
    subplot(3,1,2);
    imshow(newim);
    title('Colored Blob Mask', 'FontSize', fontSize);
	drawnow;
    
    [labeledImage, numberOfRegions] = bwlabel(newim);
	if numberOfRegions == 7
		stats = regionprops(labeledImage, 'BoundingBox', 'Centroid');
        count = count +1;
		% Delete old texts and rectangles
		if exist('hRect', 'var')
			delete(hRect);
		end
		if exist('hText', 'var')
			delete(hText);
		end
		
		% Display the original image again.
		subplot(3, 1, 3); % Switch to original image.
		hImage=subplot(3, 1, 3);
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
            thisCentroid = [thisCentroid(1)-x_correction y_correction-thisCentroid(2)]*pix2mm;
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
            if(count==1 && r==4)
                centroid4(count,:) = thisCentroid;
                filenametowrite = 'cent4.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
             if(count==1 && r==5)
                centroid5(count,:) = thisCentroid;
                filenametowrite = 'cent5.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
             end
             if(count==1 && r==6)
                centroid6(count,:) = thisCentroid;
                filenametowrite = 'cent6.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
             end
             if(count==1 && r==7)
                centroid7(count,:) = thisCentroid;
                filenametowrite = 'cent7.txt';
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
                        
                X4 = [thisCentroid;centroid4(count-1,:)];
                d4 = pdist(X4,'euclidean');    
                
                X5 = [thisCentroid;centroid5(count-1,:)];
                d5 = pdist(X5,'euclidean');
                            
                X6 = [thisCentroid;centroid6(count-1,:)];
                d6 = pdist(X6,'euclidean');
                            
                X7 = [thisCentroid;centroid7(count-1,:)];
                d7 = pdist(X7,'euclidean');
                            
                     
        if(d1<d2 && d1<d3 && d1<d4 && d1<d5 && d1<d6 && d1<d7)
            centroid1(count,:) = thisCentroid;
                filenametowrite = 'cent1.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d2<d1 && d2<d3 && d2<d4 && d2<d5 && d2<d6 && d2<d7)
            centroid2(count,:) = thisCentroid;
                filenametowrite = 'cent2.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
       
        
        if(d3<d1 && d3<d2 && d3<d4 && d3<d5 && d3<d6 && d3<d7)
            centroid3(count,:) = thisCentroid;
                filenametowrite = 'cent3.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d4<d1 && d4<d2 && d4<d3 && d4<d5 && d4<d6 && d4<d7)
            centroid4(count,:) = thisCentroid;
                filenametowrite = 'cent4.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d5<d1 && d5<d2 && d5<d3 && d5<d4 && d5<d6 && d5<d7)
            centroid5(count,:) = thisCentroid;
                filenametowrite = 'cent5.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d6<d1 && d6<d2 && d6<d3 && d6<d4 && d6<d5 && d6<d7)
            centroid6(count,:) = thisCentroid;
                filenametowrite = 'cent6.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');c;
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d7<d1 && d7<d2 && d7<d3 && d7<d4 && d7<d5 && d7<d6)
            centroid7(count,:) = thisCentroid;
                filenametowrite = 'cent7.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
            end
            
            fclose('all');
            
			hRect(r) = rectangle('Position', thisBB, 'EdgeColor', 'r', 'LineWidth', 2);
			hSpot = plot(thisCentroid(1), thisCentroid(2), 'y+', 'MarkerSize', 8, 'LineWidth', 2);
			hText(r) = text(thisBB(1), thisBB(2)-20, strcat('X: ', num2str(round(thisCentroid(1))), '    Y: ', num2str(round(thisCentroid(2)))));
			set(hText(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 8, 'Color', 'yellow');
		end
		hold off
		drawnow;
       
        
    % Calculating and plotting the body rotation(change in orientation)
        
%     if(count~=1 && r==3)  
%         a = load('cent1.txt');
%         b = load('cent2.txt');
%         c = load('cent3.txt');
%     % Reading previous coordinate of the points
%        vi1x = a(1,1);
%        vi1y = a(1,2);
%        vi2x = b(1,1);
%        vi2y = b(1,2);
%        vi3x = c(1,1);
%        vi3y = c(1,2);
% 
%     % Reading the current coordinate of the points
%       vf1x = a(count,1);
%       vf1y = a(count,2);
%       vf2x = b(count,1);
%       vf2y = b(count,2);
%       vf3x = c(count,1);
%       vf3y = c(count,2);
%     
% 
%      % Angle between 2 vectors - initial state and final state
%         theta_1 = orientation(vi2x,vi1x,vi2y,vi1y,vf2x,vf1x,vf2y,vf1y);
%         
%         theta_2 = orientation(vi3x,vi1x,vi3y,vi1y,vf3x,vf1x,vf3y,vf1y);
% 
%         theta_3 = orientation(vi3x,vi2x,vi3y,vi2y,vf3x,vf2x,vf3y,vf2y);
% 
% 
% 
%     %  Mean of all theta's
%         theta = ((theta_1 + theta_2 + theta_3)/3);
%     
%         %Plot    
%         subplot(4,1,4);
%         plot(k,theta,'r*');
%             xlim([0 260]);
%             ylim([-80 200]);
%             drawnow 
%             hold on 
%     end
    
    end
    
    
	% See if you want to bail out
	if get(hCheckbox, 'Value')
		% Finish now checkbox is checked.
		msgbox('Done with processing.');
		return;
	end
end
msgbox('Done with processing.');

%%
% Plots of the tracked points. Also to check whether the tracking has been
% done properly.

figure(2)

        a = load('cent1.txt');
        b = load('cent2.txt');
        c = load('cent3.txt');
        d = load('cent4.txt');
        e = load('cent5.txt');
        f = load('cent6.txt');
        g = load('cent7.txt');
        
subplot(3,3,1);
plot(a(:,1));
hold on
plot(a(:,2));
legend('Point1 x','Point1 y');
subplot(3,3,2);
plot(b(:,1));
hold on
plot(b(:,2));
legend('Point2 x','Point2 y');
subplot(3,3,3);
plot(c(:,1));
hold on
plot(c(:,2));
legend('Point3 x','Point3 y');
subplot(3,3,4);
plot(d(:,1));
hold on
plot(d(:,2));
legend('Point4 x','Point4 y');
subplot(3,3,5);
plot(e(:,1));
hold on
plot(e(:,2));
legend('Point5 x','Point5 y');
subplot(3,3,6);
plot(f(:,1));
hold on
plot(f(:,2));
legend('Point6 x','Point6 y');
subplot(3,3,7);
plot(g(:,1));
hold on
plot(g(:,2));
legend('Point7 x','Point7 y');

% % Identifying the respective points/retaining the marker points

A = [a(1,1) b(1,1) c(1,1) d(1,1) e(1,1) f(1,1) g(1,1)];
MIN = min(A);
MAX = max(A);

if (MIN == a(1,1))
    
end
if (MIN == b(1,1))
end
if (MIN == c(1,1))
end
if (MIN == d(1,1))
end
if (MIN == e(1,1))
end
if (MIN == f(1,1))
end
if (MIN == g(1,1))
end



% % Unit directional vectors of sides of triangle
% vi12 = [(vi2x-vi1x) (vi2y-vi1y) 0]/norm([(vi2x-vi1x) (vi2y-vi1y) 0]);
% vi13 = [(vi3x-vi1x) (vi3y-vi1y) 0]/norm([(vi3x-vi1x) (vi3y-vi1y) 0]);
% vi23 = [(vi3x-vi2x) (vi3y-vi2y) 0]/norm([(vi3x-vi2x) (vi3y-vi2y) 0]);
% vf12 = [(vf2x-vf1x) (vf2y-vf1y) 0]/norm([(vf2x-vf1x) (vf2y-vf1y) 0]);
% vf13 = [(vf3x-vf1x) (vf3y-vf1y) 0]/norm([(vf3x-vf1x) (vf3y-vf1y) 0]);
% vf23 = [(vf3x-vf2x) (vf3y-vf2y) 0]/norm([(vf3x-vf2x) (vf3y-vf2y) 0]);
% 
% 
% 
% % Identifying the respective points/retaining the marker points
% gamma_1 = radtodeg((atan2(norm(cross(vf12,vf13)),dot(vf12,vf13))));
% gamma_2 = radtodeg((atan2(norm(cross(vf12,vf23)),dot(vf12,vf23))));
% gamma_3 = radtodeg((atan2(norm(cross(vf13,vf23)),dot(vf13,vf23))));
% 
% 
% 
% % Specifying the conditions of the sides of the triangle to retain the
% % respective points
figure(3)
imshow(thisFrame);
title('Right tendon only', 'FontSize', fontSize);
hold on
 

filenametowrite = 'cent1.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
plot(centroid1(113,1),centroid1(113,2),'y+','MarkerSize',4,'LineWidth',1);

filenametowrite = 'cent7.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
plot(centroid7(113,1),centroid7(113,2),'y+','MarkerSize',4,'LineWidth',1);
 

% Draw a straight line between 2 points
A = [centroid1(113,1) centroid7(113,1)]; 
B = [centroid1(113,2) centroid7(113,2)]; 
line(A,B,'LineWidth',2)

r_r_cord = [centroid1(113,:);centroid7(113,:)];
r_r = pdist(r_r_cord,'euclidean');

% if((85<gamma_2)&&(gamma_2<95)) 
% filenametowrite = 'cent2.txt';
% fulltextfilename = fullfile(folder,filenametowrite);
% fileid = fopen(fulltextfilename,'a+');
% fprintf(fileid,'Point B');
% plot(b(k,1),b(k,2),'y+','MarkerSize',4,'LineWidth',1);
% text(b(k,1),b(k,2),'B','Color','black','FontSize',10);
% end
% 
% 
% if((60<gamma_3)&&(gamma_3<70))
% filenametowrite = 'cent3.txt';
% fulltextfilename = fullfile(folder,filenametowrite);
% fileid = fopen(fulltextfilename,'a+');
% fprintf(fileid,'Point C');
% plot(c(k,1),c(k,2),'y+','MarkerSize',4,'LineWidth',1);
% text(c(k,1),c(k,2),'C','Color','black','FontSize',10);
% end
% 
% % Calculating the Rotation matrix
% Ai = [vi12;vi13;vi23];
% Af = [vf12;vf13;vf23];
% Rot_m = pinv(Af)*Ai;
















