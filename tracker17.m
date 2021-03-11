%%% This code is just a combined version of tracker03.m and tracker05.m
%%% This code will perform real time tracking. The video is live relayed
%%% from webcam. 
%%% Live tracking of Caitlin Inchworm
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

% Capturing the plot every frame to create a movie
vwrite = VideoWriter('gcfcaptureV2','MPEG-4');
open(vwrite);

% Correction factor/values - x axis and y axis

x_correction = 0;
y_correction = 0;

% Pixel to mm conversion

pix2mm = 1;  % This is only for this particular video. Should be updated for different videos

numberOfFrames = 5000;
% Read one frame at a time and find specified color

while true
    
    if k == 1
		set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
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
        
        % Create a file in the folder in which you are working to save the
        % datas in the particular file.
        
%         a = load('cent1.txt');
%         b = load('cent2.txt');
%         c = load('cent3.txt');
%         d = load('cent4.txt');
%         e = load('cent5.txt');
%         f = load('cent6.txt');
%         g = load('cent7.txt');  
	end
    
    
    thisFrame=cam.snapshot;
    thisFrame = imresize(thisFrame,[270,480]); % Changing the resolution
    newim = createMaskInchworm2(thisFrame);
    % Filter out small blobs
    newim = bwareaopen(newim, 10);
    % Fill in holes
    newim = imfill(newim, 'holes');
    hImage=subplot(3, 1, 1);
    % Display it.
	imshow(thisFrame);
	axis on;
	caption = sprintf('Original image, frame #%d 0f %d', k, numberOfFrames);
	title(caption, 'FontSize', fontSize);
	drawnow;
%     subplot(3,1,2);
%     imshow(newim);
%     title('Colored Blob Mask', 'FontSize', fontSize);
% 	drawnow;
    
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
		
		% Displaying original image again.
		subplot(3, 1, 3); 
		hImage=subplot(3, 1, 3);
		imshow(thisFrame);
		axis on;
		hold on;
		caption = sprintf('%d blobs found in frame #%d 0f %d', numberOfRegions, k, numberOfFrames);
		title(caption, 'FontSize', fontSize);
		drawnow;
%       captureFrame = getframe(gcf);
%       writeVideo(vid,captureFrame);
        %This is a loop to bound the colored objects in a rectangular box.
		for r = 1 : numberOfRegions
			% Find location for this blob.
			thisBB = stats(r).BoundingBox;
			thisCentroid = stats(r).Centroid;
            thisCentroid = [thisCentroid(1)-x_correction thisCentroid(2)]*pix2mm;
            
            % Nearest neighbour
            if(count==1 && r==1)
                centroid1(count,:) = thisCentroid;
%                 filenametowrite = 'cent1.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
            if(count==1 && r==2)
                centroid2(count,:) = thisCentroid;
%                  filenametowrite = 'cent2.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
            if(count==1 && r==3)
                centroid3(count,:) = thisCentroid;
%                 filenametowrite = 'cent3.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
            if(count==1 && r==4)
                centroid4(count,:) = thisCentroid;
%                 filenametowrite = 'cent4.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
            end
             if(count==1 && r==5)
                centroid5(count,:) = thisCentroid;
%                 filenametowrite = 'cent5.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
             end
             if(count==1 && r==6)
                centroid6(count,:) = thisCentroid;
%                 filenametowrite = 'cent6.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
             end
             if(count==1 && r==7)
                centroid7(count,:) = thisCentroid;
%                 filenametowrite = 'cent7.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
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
%                 filenametowrite = 'cent1.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d2<d1 && d2<d3 && d2<d4 && d2<d5 && d2<d6 && d2<d7)
            centroid2(count,:) = thisCentroid;
%                 filenametowrite = 'cent2.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
       
        
        if(d3<d1 && d3<d2 && d3<d4 && d3<d5 && d3<d6 && d3<d7)
            centroid3(count,:) = thisCentroid;
%                 filenametowrite = 'cent3.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d4<d1 && d4<d2 && d4<d3 && d4<d5 && d4<d6 && d4<d7)
            centroid4(count,:) = thisCentroid;
%                 filenametowrite = 'cent4.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d5<d1 && d5<d2 && d5<d3 && d5<d4 && d5<d6 && d5<d7)
            centroid5(count,:) = thisCentroid;
%                 filenametowrite = 'cent5.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d6<d1 && d6<d2 && d6<d3 && d6<d4 && d6<d5 && d6<d7)
            centroid6(count,:) = thisCentroid;
%                 filenametowrite = 'cent6.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');c;
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
        
        if(d7<d1 && d7<d2 && d7<d3 && d7<d4 && d7<d5 && d7<d6)
            centroid7(count,:) = thisCentroid;
%                 filenametowrite = 'cent7.txt';
%                 fulltextfilename = fullfile(folder,filenametowrite);
%                 fileid = fopen(fulltextfilename,'a+');
%                 fprintf(fileid,'%d   %d \n',thisCentroid(1),thisCentroid(2));
        end
            end
            
%             fclose('all');
            
			hRect(r) = rectangle('Position', thisBB, 'EdgeColor', 'r', 'LineWidth', 2);
			hSpot = plot(thisCentroid(1), thisCentroid(2), 'y+', 'MarkerSize', 7, 'LineWidth', 2);
			hText(r) = text(thisBB(1), thisBB(2)-20, strcat('X: ', num2str(round(thisCentroid(1))), '    Y: ', num2str(round(thisCentroid(2)))));
			set(hText(r), 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 8, 'Color', 'yellow');
%             captureFrame = getframe(gcf);
%             writeVideo(vid,captureFrame);
		end
		hold off
		drawnow;
    frame = getframe(gcf);
    writeVideo(vwrite,frame);
    end
    
    
	% To stop the process in between.
	if get(hCheckbox, 'Value')
		msgbox('Quit processing');
		return;
	end
end
msgbox('Quit processing');
close(vwrite);

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


figure(3)
imshow(thisFrame);
title('Retaining marker points - Both actuated', 'FontSize', fontSize);
hold on
 

filenametowrite = 'cent1.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
plot(centroid1(158,1),centroid1(158,2),'y+','MarkerSize',4,'LineWidth',1);


filenametowrite = 'cent7.txt';
fulltextfilename = fullfile(folder,filenametowrite);
fileid = fopen(fulltextfilename,'a+');
plot(centroid7(158,1),centroid7(158,2),'y+','MarkerSize',4,'LineWidth',1);


% Draw a straight line between 2 points
A = [centroid1(158,1) centroid7(158,1)]; 
B = [centroid1(158,2) centroid7(158,2)]; 
line(A,B,'LineWidth',2)

sigma_min_cord = [centroid1(158,:);centroid7(158,:)];
sigma_min = pdist(sigma_min_cord ,'euclidean');












