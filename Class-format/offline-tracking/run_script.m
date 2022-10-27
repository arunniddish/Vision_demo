clear all;
clc;

% Adding dependencies
addpath('../');

% Video file information
input_vid_filename = 'Gait F/Trial_1.mp4';
output_vid_filename = 'Gait F/Trial_1_gcf.mp4';

% Instantiate a video objects for this video.
params.vread = VideoReader(input_vid_filename);
params.vwrite = VideoWriter(output_vid_filename,'MPEG-4');
open(params.vwrite);

% Tracking parameters
params.number_of_markers = 4;

% Initilaize
tracker_obj = OfflineTracking(params);

% Call tracking function 
output_data =  tracker_obj.tracking();