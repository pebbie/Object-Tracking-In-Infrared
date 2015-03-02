%GetGTPlots:
%This function extracts groundtruth data or extended groundtruth data.
%Normal groundtruth data accompanies that OCTBVS database, Dataset 1 and is
%sufficient for Single_Track. However, extended ground_truth data is needed
%for NoDetect_Track. Extended ground_truth data involved manually tacking
%each object detection in a frame as new or old. GetGTPlots extracts the
%ground truth location data for use by other functions. Note: Credit for
%basis of the ground truth extraction code goes to the OCTBVS group
%(cited in the accompanying paper). That code snippet was modified and
%expanded to create GetGTPlots. It was also rebuilt to accomodate the
%modified extended ground_truth data files. This version of GetGTPlots
%extracts extended ground truth data.

function [GTPlot] = GetGTPlots(~)
GTPlot = [];

%cd('C:\Users\Pushpak\Documents\MATLAB\OTCBVS\Dataset 1');
%First access the Dataset 1/5 Directory list for the selected folder:
%cd(c_folder);

%Grab the Groundtruth file with associated sequence
fid = fopen('groundTruth_Enhanced.txt');
    
%Cycle through comments in the file:
line = fgets(fid);
while line(1) == '%'
    line = fgets(fid);
end
    
    
%Read the number of images in the video
numImages = sscanf(line, '%d', 1);
start = 1;
    
for i=1:numImages
        
    %Load the image name and number of boxes
    imageName = fscanf(fid, '%c',13);
    numBoxes = fscanf(fid, '%d', 1);

    %Display the image
    fname = fullfile(imageName);
    Im = imread(fname);

    %Load the ground truth boxes
    for j=1:numBoxes
        tmp = fscanf(fid, '%c',1); %% [space](
        Object_Status = fscanf(fid, '%c',1);
        tmp = fscanf(fid, '%c',1); %% [space](
        coords = fscanf(fid, '%d %d %d %d');
        tmp = fscanf(fid, '%c',1); %% )
        ulX=coords(1); ulY=coords(2);
        lrX=coords(3); lrY=coords(4);
        GTPlot{i,j}.c = [ulX ulY abs(ulX-lrX) abs(ulY-lrY)];
        GTPlot{i,j}.x = [ulX lrX lrX ulX ulX];
        GTPlot{i,j}.y = [ulY ulY lrY lrY ulY];
    end

    tmp = fgetl(fid); %% get until end of line

end

    fclose(fid);

end





