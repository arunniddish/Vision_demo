%%% This code is same as tracker 07
%%% This code is to track Caitlin's video - Inch Worm
%%% The Inch worm has 7 markers on it. Therefore 7 different markers are
%%% tracked seperately. 
%%% Natural length


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
baseFileName = 'InchWorm_Trim.mp4';  % Video to track
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



% Read one frame at a time and find specified color
for k = 1:numberOfFrames
    
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
    newim = createMaskInchWorm(thisFrame);
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
            
            %Identify RED color LED 
    
    newimRED = createMaskRED(thisFrame);
    % Filter out small blobs
    newimRED = bwareaopen(newimRED, 10);
    % Fill in holes
    newimRED = imfill(newimRED, 'holes');
    [labeledImage_RED, numberOfRegions_RED] = bwlabel(newimRED);
    if exist('numberOfRegions_1','var')
        stats_RED = regionprops(labeledImage_RED, 'BoundingBox', 'Centroid');
		% Delete old texts and rectangles
		if exist('hRect_RED', 'var')
			delete(hRect_RED);
		end
		if exist('hText_RED', 'var')
			delete(hText_RED);
        end
   
        %This is a loop to bound the colored objects in a rectangular box.
		for r = 1 : numberOfRegions_RED
			% Find location for this blob.
			thisBB_RED = stats_RED(r).BoundingBox;
			thisCentroid_RED = stats_RED(r).Centroid;
        
        hRect_RED(r) = rectangle('Position', thisBB_RED, 'EdgeColor', 'r', 'LineWidth', 2);
		hSpot_RED = plot(thisCentroid_RED(1), thisCentroid_RED(2), 'y+', 'MarkerSize', 8, 'LineWidth', 2);
        str_RED = 'Red light is activated';
 		hText_RED(r) = text(thisBB_RED(1), thisBB_RED(2)-20, str_RED);
		set(hText_RED(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 8, 'Color', 'yellow');
        end
    end

    %Identify BLUE color LED 
    
    newimBLUE = createMaskRED(thisFrame);
    % Filter out small blobs
    newimBLUE = bwareaopen(newimBLUE, 10);
    % Fill in holes
    newimBLUE = imfill(newimBLUE, 'holes');
    [labeledImage_BLUE, numberOfRegions_BLUE] = bwlabel(newimBLUE);
    if exist('numberOfRegions_BLUE','var')
        stats_BLUE = regionprops(labeledImage_BLUE, 'BoundingBox', 'Centroid');
		% Delete old texts and rectangles
		if exist('hRect_BLUE', 'var')
			delete(hRect_BLUE);
		end
		if exist('hText_BLUE', 'var')
			delete(hText_BLUE);
        end
   
        %This is a loop to bound the colored objects in a rectangular box.
		for r = 1 : numberOfRegions_BLUE
			% Find location for this blob.
			thisBB_BLUE = stats_BLUE(r).BoundingBox;
			thisCentroid_BLUE = stats_BLUE(r).Centroid;
        
        hRect_BLUE(r) = rectangle('Position', thisBB_BLUE, 'EdgeColor', 'r', 'LineWidth', 2);
		hSpot_BLUE = plot(thisCentroid_BLUE(1), thisCentroid_BLUE(2), 'y+', 'MarkerSize', 8, 'LineWidth', 2);
        str_BLUE = 'Blue light is activated';
 		hText_BLUE(r) = text(thisBB_BLUE(1), thisBB_BLUE(2)-20, str_BLUE);
		set(hText_BLUE(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 8, 'Color', 'yellow');
        end
    end
    
    %Identify GREEN color LED 
    
    newimGREEN = createMaskRED(thisFrame);
    % Filter out small blobs
    newimGREEN = bwareaopen(newimGREEN, 10);
    % Fill in holes
    newimGREEN = imfill(newimGREEN, 'holes');
    [labeledImage_GREEN, numberOfRegions_GREEN] = bwlabel(newimGREEN);
    if exist('numberOfRegions_1','var')
        stats_GREEN = regionprops(labeledImage_GREEN, 'BoundingBox', 'Centroid');
		% Delete old texts and rectangles
		if exist('hRect_GREEN', 'var')
			delete(hRect_GREEN);
		end
		if exist('hText_GREEN', 'var')
			delete(hText_GREEN);
        end
   
        %This is a loop to bound the colored objects in a rectangular box.
		for r = 1 : numberOfRegions_GREEN
			% Find location for this blob.
			thisBB_GREEN = stats_GREEN(r).BoundingBox;
			thisCentroid_GREEN = stats_GREEN(r).Centroid;
        
        hRect_GREEN(r) = rectangle('Position', thisBB_GREEN, 'EdgeColor', 'r', 'LineWidth', 2);
		hSpot_GREEN = plot(thisCentroid_GREEN(1), thisCentroid_GREEN(2), 'y+', 'MarkerSize', 8, 'LineWidth', 2);
        str_GREEN = 'Red light is activated';
 		hText_GREEN(r) = text(thisBB_GREEN(1), thisBB_GREEN(2)-20, str_GREEN);
		set(hText_GREEN(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 8, 'Color', 'yellow');
        end
    end
            
		end
		hold off
		drawnow;
        
    
        
        
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
title('Retaining marker points', 'FontSize', fontSize);
hold on
 

filenametowrite = 'cent1.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
plot(a(1,1),a(1,2),'y+','MarkerSize',4,'LineWidth',1);

filenametowrite = 'cent7.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
plot(g(1,1),g(1,2),'y+','MarkerSize',4,'LineWidth',1);
 

% Draw a straight line between 2 points
A = [a(1,1) g(1,1)]; 
B = [a(1,2) g(1,2)]; 
line(A,B,'LineWidth',2)

length_cord = [a(1,:);g(1,:)];
n_len = pdist(length_cord,'euclidean');

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
















