%Handoffs:
%This function extracts extended groundtruth data.Extended ground_truth 
%data is needed for NoDetect_Track. This function extracts the manually
%input labels in the extended ground truth data. The labels are returned in
%a mask format that corresponds to the ground truth matrices produced by
%GetGTPlots and GetGTPlots1

function [handoffs] = Handoff(~)

handoffs = [];
%First access the Dataset 1/5 Directory list for the selected folder:
%cd('C:\Users\Pushpak\Documents\MATLAB\OTCBVS\Dataset 1');
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
    
for i=1:numImages
        
    %Load the image name and number of boxes
    imageName = fscanf(fid, '%c',13);
    numBoxes = fscanf(fid, '%d', 1);

    fname = fullfile(imageName);

    for j=1:numBoxes
        tmp = fscanf(fid, '%c',1); %% [space](
        Object_Status = fscanf(fid, '%c',1);
        
        if Object_Status == 'N'
            handoffs(i,j) = 1;
        else
            handoffs(i,j) = 0;
        end
        
        tmp = fscanf(fid, '%c',1); %% [space](
        tmp = fscanf(fid, '%d %d %d %d');
        tmp = fscanf(fid, '%c',1); %% )

    end

    tmp = fgetl(fid); %% get until end of line

end

    fclose(fid);

end

