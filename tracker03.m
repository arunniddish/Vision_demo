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
%         centroid1 = zeros(1000,2);
%         centroid2 = zeros(1000,2);
%         centroid3 = zeros(1000,2);
	end
    
    
    thisFrame=cam.snapshot;
    thisFrame = imresize(thisFrame,[270,480]); % Resolution of image is 270 x 480 pixel
    newim = createMask(thisFrame);
    % Filter out small blobs
    newim = bwareaopen(newim, 200);
    % Fill in holes
    newim = imfill(newim, 'holes');
    hImage=subplot(3, 1, 1);
    % Display it.
	imshow(thisFrame);
	axis on;
	caption = sprintf('Original RGB image, frame #%d 0f %d', k, k);
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
		caption = sprintf('%d blobs found in frame #%d 0f %d', numberOfRegions, k, k);
		title(caption, 'FontSize', fontSize);
		drawnow;
        
        %This is a loop to bound the colored objects in a rectangular box.
		for r = 1 : numberOfRegions
			% Find location for this blob.
			thisBB = stats(r).BoundingBox;
			thisCentroid = stats(r).Centroid;
       if(numberOfRegions==7)
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











