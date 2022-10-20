clear all;
clc;

addpath('../');

input_vid_filename = 'Gait F/Trial_1.mp4';
output_vid_filename = 'Gait F/Trial_1_gcf.mp4';

% Instantiate a video objects for this video.
params.vread = VideoReader(input_vid_filename);
params.vwrite = VideoWriter(output_vid_filename,'MPEG-4');

% Tracking parameters
params.number_of_markers = 4;

tracker_obj = OfflineTracking(params);
output_data =  tracker_obj.tracking();

