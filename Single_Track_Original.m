%-------------------------------------------------------------------------%
%Single_Track  - Generalized Matlab Script capable of implementing
%single color (in this case IR) tracking techniques. Currently
%contains image processing to remove noise and to emphasize detected
%objects. This includes background averaging, and subsequent subtraction
%from the current frame to detect foreground elements.
%
%Technique 1: Foreground objects are maintained in a tracking list. Object 
%uncertainty starts as moderate. If the object is not identified in
%subsequent frames, the uncertainty increases. If the object is identified
%the uncertainty decreases. Objects that are stably tracked over time
%therefore have lower uncertainty and objects that are not will have
%uncertainty. If uncertainty passes a certain threshold the object is
%removed from the track list. If an object on the track list cannot be
%associated but has low uncertainty, it is maintained and the track box is
%placed at the predicted object location. 
%
%Technique 2: Implements a MEANSHIFT technique to find distribution centers 
%identified using technique 1. This can be turned on or off. The
%MEANSHIFT identification is done in the seperate Matlab function camshift.m.
%Camshift1 or 2 can be used to track the distribution on the original image,
%or more processed images at any stage of the alogirthm. This provided
%antoher 'lever' to pull to affect tracking performance. It would be
%possible to use the MEANSHIFT based algorithms to track the objects on its
%own but this was less effective than technique 1. These identified points
%are included as a demonstration of distributional centers identified using
%technique 1.
%-------------------------------------------------------------------------%

%addpath('C:\Users\aa\Downloads\APOGEE\Pushpak-Object Tracking\Pushpak-Object tracking-code\Attachments_Complete');
%cd('C:\Users\Pushpak\Documents\MATLAB\OTCBVS\Dataset 1');

%Select the desired folder
c_folder = '00001';

%First access the Dataset 1/5 Directory list for the selected folder:
cd(c_folder);
%clear all;
close all;
disp('Program Start...');

%Variables:
d_list = dir('*.bmp');  %Directory listing of all frames
background_status = 1;  %Sets background status (0 for bg, 1 for fg)
bgf_count = 1;          %Background frame count  
ms_win = zeros(1,4);    %Meanshift search window - Starts empty
ms_x = 0;               %Meanshift x plot location
ms_y = 0;               %Meanshift y plot location
t_list = [];            %Track list keeps record of each track, as well as 
                        %the previous tracked location and a rough
                        %prediction of where it will be with the following 
                        %formatting: 
                        %[C_x C_y BB_x BB_y Prev_x Prev_y Pre_x Pre_y]
b_pix = 15;             %Measure of acceptable closeness to the boundary
uncertainty = 33;       %Measure of uncertainty of a track. Increases if 
                        %is not detected in frames and decreases to 0 if
                        %the object is redetected
o_size = 25;            %Minimum detectable object size
cs_list = [];           %CAM-Shift Tracking List

%Background value: Updates with all non-object frames. Assumes first frame 
%is background
I_p = double(imread(d_list(2).name)) - double(imread(d_list(1).name));

for n1 = 1:size(I_p,1)
    for n2 = 1:size(I_p,2)
        if I_p(n1,n2) > 0
            I_p(n1,n2) = 0;
        end
    end
end
I_bg = double(imread(d_list(1).name)) + I_p; 
                              
disp(['Expected Playback Time: ', num2str(numel(d_list)*.4/60),' mins']);


%Main Body: Plays videos and performs tracking using various methods

tic
for i = 2: numel(d_list)
    
    I_curr = double(imread(d_list(i).name));
    %Previous can be set as the previous frame for motion templating I_prev
    %= imread(d_list(i).name);
    
%-------------------------------------------------------------------------%
%Preprocessing: 
%Performs background subtraction, binarization off of an empirically
%defined threshold, erosion to remove small noise objects, and finally
%displays the processed image as well as the original image in seperate
%figures
    
    %As a default, assume objects are not detected (background frame)
    background_status = 1;
    
    %Take the difference between frame and bg to find foreground
    I_diff = abs(I_curr-I_bg);      
    %Binarize with an experimentally determined high threshold
    
    I_bdiff = I_diff > 65;
    %I_bdiff = I_diff > 30;
    %I_bdiff = I_diff>50;
    
    %Erode to remove noise
    erodesize = 1;
    I_db = imerode(I_bdiff,ones(erodesize,erodesize));
    
    %Display the image frame and processing frame
    figure(1),imshow(d_list(i).name); 
    hold on;

%-------------------------------------------------------------------------%
%Creating the track list:
%If an object is detected then we dilate the remaining objects to increase
%the likelihood of internally connecting discrete objects. Image labeling
%is used to identify all seperate objects, and bounding boxes are
%determined for each object. If a list already exists, clear the previous
%location elements and replace with a '-1' marker. 

    %If we have any objects:
    if max(max(I_db))>0
        
        %If objects, then this is not background
        background_status = 0;
        %Dilate to keep binarized image as internally connected
        
        dy = 12;
        dx = 8;
        I_db = imdilate(I_db,ones(dy,dx));
        %I_db = bwmorph(I_db, 'majority');
        
        figure(2),imshow(I_db,[0 1])
        hold on;
        
        %Label the image regions 
        [I_l num_l] = bwlabel(I_db,8);
        L_stats = regionprops(I_l,'BoundingBox');
        
        %Uses Region Labeling Approach, creating a bounding box around
        %sections. Plots the bounding box around each region
        
        diffs = [];  %Assigns best guess for detections
        
        if isempty(t_list)==0
            
            %First, remove stored 'previous points' in the t_list
            for n = 1:size(t_list,1)
                t_list(n,5) = -1;
                t_list(n,6) = -1;
            end
            
            %Then, create an ordered difference list between current points
            %and previously stored track points
            for n = 1:num_l
                if (L_stats(n).BoundingBox(3)+L_stats(n).BoundingBox(4)) > o_size
                    
                    ix = L_stats(n).BoundingBox(1);
                    iy = L_stats(n).BoundingBox(2);
                    
                    %Check if this object is on the track list by creating a
                    %difference list between current locations and all
                    %predicted locations in the track list:
                    
                    diff_list = abs(ix - t_list(:,7)) + abs(iy - t_list(:,8));
                    

                else
                    diff_list = ones(size(t_list,1),1);
                end
                
                diffs = [diffs diff_list];
            end                    
            
            %The number of unassociated track points. Initally the length
            %of the list
            num_uatp = size(t_list,1);
            used_list = zeros(size(diffs,2),1);
            
            %Next, associate detected points with old track points. If not
            %associated points are found, then add new track points.
 
            for n = 1:size(diffs,2)  
                %Find the closest associated predicted track to a new
                %detection
                %Reset idx1:
                idx1 = [];
                
                if (size(diffs,1) > 1 && size(diffs,2) > 1)
                    [m1 idx1] = min(diffs);
                    [m2 idx2] = min(m1);
                else
                    if size(diffs,1) == 1
                        [m1 idx2] = min(diffs);
                        idx1(idx2) = 1;
                    end
                    if size(diffs,2) == 1
                        idx2 = 1;
                        [m1 idx1(idx2)] = min(diffs);
                    end    
                end
                
                y1 = idx1(idx2);
                x1 = idx2;
                
                %If the points still are fairly well associated and also
                %there are track points left to associate. If a point has
                %been 'missing' decrease the threshold of detection for
                %surrounding points by increasing the uncertainty marker
                if num_uatp > 0 && diffs(y1,x1) < ...
                   (120 + t_list(y1,9)) 
                    %Update the current track point and the previous point
                    t_list(y1,5) = t_list(y1,1); 
                    t_list(y1,6) = t_list(y1,2);
                    t_list(y1,1) = L_stats(x1).BoundingBox(1);
                    t_list(y1,2) = L_stats(x1).BoundingBox(2);
                    t_list(y1,3) = L_stats(x1).BoundingBox(3); 
                    t_list(y1,4) = L_stats(x1).BoundingBox(4);
                    %Decrease the uncertainty once a detection:
                    t_list(y1,9) = max([(t_list(y1,9) - uncertainty) 0]);
                    
                    if isempty(diffs) == 0
                        %Fill the used row
                        diffs(y1,:) = 1000 * ones(1,size(diffs,2));
                        %Fill the used column
                        diffs(:,x1) = 1000 * ones(size(diffs,1),1);
                    end
                    
                    %An old track point has now been associated:
                    num_uatp = num_uatp - 1; 
                    used_list(x1) = 1;
                end
            end
            
            %If we are out of track points to associate or no close
            %track points are available, then we create a new track
            %point based on the used_list
            for n = 1:size(diffs,2)  
                if used_list(n) == 0
                    ix = L_stats(n).BoundingBox(1);
                    iy = L_stats(n).BoundingBox(2);
                    sx = L_stats(n).BoundingBox(3);
                    sy = L_stats(n).BoundingBox(4);
                    
                    mb = [abs(ix - 1) abs(iy - sy - 1)...
                          abs(ix + sx -size(I_db,2)) abs(iy - size(I_db,1))];
                    [tm pos] = min(mb);
                    
                    s_tlist = size(t_list,1);
                    %Based on which side wall is closest, fill the t_list
                    %These points have moderate uncertainty (50)
                        if pos == 1
                            t_list(s_tlist+1,:) = [ix iy sx sy 1 iy 0 0 50];
                        end
                        if pos == 2
                            t_list(s_tlist+1,:) = [ix iy sx sy ix 1 0 0 50];
                        end
                        if pos == 3
                            t_list(s_tlist+1,:) = [ix iy sx sy ...
                                                   size(I_db,2) iy 0 0 50];
                        end
                        if pos == 4
                            t_list(s_tlist+1,:) = [ix iy sx sy ix ...
                                                   size(I_db,1) 0 0 50];
                        end
                end
            end
                
            %Next, we delete all track points that are no longer detected
            %and are expected to be gone
            %Uses b_pix and uncertainty from the Variable List
                            
            for n = size(t_list,1):-1:1
                if t_list(n,5) == -1 && t_list(n,6) == -1
                    
                    %If close to the boundary and missing, delete:
                    %Alternatively, if the track has not been found in too
                    %long, then the uncertainty is too high and it will be
                    %deleted:
                    if t_list(n,7) < b_pix ||...
                       t_list(n,8) - t_list(n,4) < b_pix ||...
                       t_list(n,7) + t_list(n,3) > (size(I_db,2)-b_pix) ||...
                       t_list(n,8) > (size(I_db,1)-b_pix) ||...
                       t_list(n,9) > 50
                        
                       t_list(n,:) = [];
                    else
                    %Else guess its position based on a simple prediction
                    %and hope it comes back. Increase the uncertainty,
                    %assume the bounding box stays the same, make the
                    %current location the last predicted location, and the
                    %previous location the current location
                       t_list(n,5) = t_list(n,1);
                       t_list(n,6) = t_list(n,2);
                       t_list(n,1) = t_list(n,7);
                       t_list(n,2) = t_list(n,8);
                       t_list(n,9) = t_list(n,9) + uncertainty;                                       
                    end
                end
            end
            
            %If an object has become too small... then delete it
            for n = size(t_list,1):-1:1
                if (t_list(n,3) + t_list(n,4)) < (o_size+2)
                    t_list(n,:) = [];
                end
            end
        else
            %If the Track list is empty. Then add all tracks from scratch!
            for n = 1:num_l
                if (L_stats(n).BoundingBox(3)+L_stats(n).BoundingBox(4)) > o_size
                    ix = L_stats(n).BoundingBox(1);
                    iy = L_stats(n).BoundingBox(2);
                    sx = L_stats(n).BoundingBox(3);
                    sy = L_stats(n).BoundingBox(4);
                    
                    mb = [abs(ix - 1) abs(iy - sy - 1)...
                          abs(ix + sx -size(I_db,2)) abs(iy - size(I_db,1))];
                    [tm pos] = min(mb);
                    
                    %Based on which side wall is closest, fill the t_list
                    %These points are high uncertainty (100) - Newly 
                    %created points become more certain as they have been 
                    %visible for a longer time. When a list of all new
                    %track points is created at once, it is particularly
                    %uncertain
                        if pos == 1
                            t_list(n,:) = [ix iy sx sy 1 iy 0 0 100];
                        end
                        if pos == 2
                            t_list(n,:) = [ix iy sx sy ix 1 0 0 100];
                        end
                        if pos == 3
                            t_list(n,:) = [ix iy sx sy size(I_db,2) iy 0 0 100];
                        end
                        if pos == 4
                            t_list(n,:) = [ix iy sx sy ix size(I_db,1) 0 0 100];
                        end
                end
            end
            cs_list = t_list;
        end
%-------------------------------------------------------------------------%
%Tracklist Updating:  
%Now that the tracklist has been populated, update it!
    if isempty(t_list) == 0
        t_list = tlupdate(t_list,I_diff);
    end
%-------------------------------------------------------------------------%
%Display Tracks Section:
%        for n = 1:num_l
        disp(['SIZE: ', num2str(size(t_list,1))]);
        for n=1:size(t_list,1)
                     
                ix = t_list(n,1);
                iy = t_list(n,2);
                sx = t_list(n,3);
                sy = t_list(n,4);
                                
                %Corresponds to Technique 1:
                if t_list(n,9) < uncertainty  
                    figure(1), plot([ix ix ix+sx ix+sx ix],...
                                    [iy iy+sy iy+sy iy iy], 'r');
                    hold on;
                else
                    figure(1), plot([ix ix ix+sx ix+sx ix],...
                                     [iy iy+sy iy+sy iy iy], 'y');
                    hold on;
                end
        end
        
        ms_list = t_list;
        
        for n = 1:size(t_list,1)   
                %Corresponds to Technique 2:
                
                    ix = t_list(n,1);
                    iy = t_list(n,2);
                    sx = t_list(n,3);
                    sy = t_list(n,4);
                    
                    %initialize a start window:
                    ms_win(n,1) = round(ix);
                    ms_win(n,2) = round(iy);
                    ms_win(n,3) = round(sx);
                    ms_win(n,4) = round(sy);
                    
                    %Apply camshift tracking - gets centroids and new
                    %window for next iteration
                    [ms_win(n,:) ms_x ms_y] = camshift(I_db,ms_win(n,:));
                    
                    ix = ms_win(n,1);
                    iy = ms_win(n,2);
                    sx = ms_win(n,3);
                    sy = ms_win(n,4);
                    
                    figure(1), plot([ix ix ix+sx ix+sx ix],...
                                    [iy iy+sy iy+sy iy iy], 'b');
                %else
                    %Apply camshift tracking - gets centroids and new
                    %window for next iteration
                %    [ms_win(n,:) ms_x ms_y] = camshift(I_db,ms_win(n,:));   
                %end
                %figure(1), plot(ms_x,ms_y, 'bx','MarkerSize',20);
                hold on;
        end
    end
    
%-------------------------------------------------------------------------%
%Background Updating:
%Improve background estimate over time by averaging in new background
%estimates
    %Update the background with curr_frame if a background frame
    if background_status == 1
        bgf_count = bgf_count + 1;
        I_bg = bgupdate(I_bg,I_curr,bgf_count);
    else
    %If not a background frame, update with our best guess at bg
        if i >=2
            I_p = double(imread(d_list(i).name)) - double(imread(d_list(i-1).name));

                for n1 = 1:size(I_p,1)
                    for n2 = 1:size(I_p,2)
                        if I_p(n1,n2) > 0
                            I_p(n1,n2) = 0;
                        end
                    end
                end
            I_nbg = double(imread(d_list(i-1).name)) + I_p; 
            bgf_count = bgf_count + 1;
            I_bg = bgupdate(I_bg,I_nbg,bgf_count);
        end
    end

    %Pause to allow the user to view the frame:
    pause(1);

end
toc

%Return to the Dataset 1/5 Directory
cd ..;
