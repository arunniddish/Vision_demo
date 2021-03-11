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
baseFileName = 'Magenta.mp4';  % Video to track
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
	end
    
    
    thisFrame=read(videoObject,k);
    newim = createMask(thisFrame);
    % Filter out small blobs
    newim = bwareaopen(newim, 200);
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
            if(k==1 && r==1)
                centroid1(k,:) = thisCentroid;
            end
            if(k==1 && r==2)
                centroid2(k,:) = thisCentroid;
            end
            if(k~=1)
                X1 = [thisCentroid;centroid1(k-1,:)];
                d1 = pdist(X1,'euclidean');
                
                X2 = [thisCentroid;centroid2(k-1,:)];
                d2 = pdist(X2,'euclidean');
            
                
           
        if(d1<d2)
            centroid1(k,:) = thisCentroid;
                filenametowrite = 'cent1.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d1>d2)
            centroid2(k,:) = thisCentroid;
                filenametowrite = 'cent2.txt';
                fulltextfilename = fullfile(folder,filenametowrite);
                fileid = fopen(fulltextfilename,'a+');
                fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
            end
                
                
            % Writing the coordinates of centroid to a file
%             if r == 1
%                 filenametowrite = 'cent1.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
%             end
%             if r == 2
%                 filenametowrite = 'cent2.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
%             end
            
			hRect(r) = rectangle('Position', thisBB, 'EdgeColor', 'r', 'LineWidth', 2);
			hSpot = plot(thisCentroid(1), thisCentroid(2), 'y+', 'MarkerSize', 10, 'LineWidth', 2)
			hText(r) = text(thisBB(1), thisBB(2)-20, strcat('X: ', num2str(round(thisCentroid(1))), '    Y: ', num2str(round(thisCentroid(2)))));
			set(hText(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
		end
		hold off
		drawnow;
    end
    
    
      	% See if they want to bail out
	if get(hCheckbox, 'Value')
		% Finish now checkbox is checked.
		msgbox('Done with demo.');
		return;
	end
end
msgbox('Done with demo.');

a = load('cent1.txt');
b = load('cent2.txt');
figure(2)
plot(a(:,1));
hold on
plot(a(:,2));
hold on
plot(b(:,1));
hold on
plot(b(:,2));
legend('Point1 x','Point1 y','Point2 x','Point2 y');













