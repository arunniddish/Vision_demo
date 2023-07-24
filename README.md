# Vision
Documentation:
https://bama365-my.sharepoint.com/:w:/g/personal/anmahendran_crimson_ua_edu/EbZkQT3FR9FHoAhuxtGwDL8BPs3d6qBJrz4_uiu66-wsTw?email=anmahendran%40crimson.ua.edu&e=4%3A8iQhp5&at=9&CID=780bc6bb-2f8f-cdf1-0515-560172756d4c

All .mat files consists of tracking data.\
Capture rate at 30fps.

The files are 'n x 48' matrix.\
n --> Number of frames.

Matlab commands used to structure the data:

theta_G(k,:) = reshape(Rot_G,[1,9]);   % Rotation matrix w.r.t 1st frame\
trans_G(k,:) = T_G';                   % Translation matrix w.r.t 1st frame

% Collectively writing all the tracked points in a single variable.\
all_pt = cat(2,centroid1,centroid2,centroid3,centroid4,centroid5,centroid6,...\
             centroid7,centroid8,theta,trans,theta_G,trans_G);
             
They are structured as follows:

Columns             Respective Datas\
1  to 3  (1x3)      Marker 1 [x,y,z]\
4  to 6  (1x3)      Marker 2 [x,y,z]\
7  to 9  (1x3)      Marker 3 [x,y,z]\
10 to 12 (1x3)      Marker 4 [x,y,z]\
13 to 15 (1x3)      Marker 5 [x,y,z]\
16 to 18 (1x3)      Marker 6 [x,y,z]\
19 to 21 (1x3)      Marker 7 [x,y,z]\
22 to 24 (1x3)      Marker 8 [x,y,z]\
25 to 33 (1x9)      Rotation matrix Intermediate frames [Reshaped from 3x3 to 1x9]\
34 to 36 (1x3)      Translation Matrix Intermediate frames [[x,y,z]Reshaped from 3x1 to 1x3] \
37 to 45 (1x9)      Rotation matrix Global frame [Reshaped from 3x3 to 1x9]\
46 to 48 (1x3)      Translation Matrix Global frame [[x,y,z]Reshaped from 3x1 to 1x3]

Pixel to mm conversion: \
1 pixel(x-axis) = 2.3921 mm  (630x400)\
1 pixel(y-axis) = 2.3874 mm  (630x400)





