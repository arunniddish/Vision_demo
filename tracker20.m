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
%marker = zeros(numberOfFrames,2);

% Initializing serial port
device = serialport('COM10',9600);

% Initializing camera
cam = webcam;

thisFrame = cam.snapshot;
thisFrame = imresize(thisFrame,[270,480]); % Changing the resolution
newim = createMaskpoor(thisFrame);

  




