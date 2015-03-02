%-------------------------------------------------------------------------%
%NoDetect_Track: Identifies object positions initally using groundtruth
%data. After this tracking approaches are employed. This method is in 
%contrast to Single_Track which uses image preprocessing toidentify targets
%
%Methodology: This method provides a clean evaluation technique of tracker 
%performance fully divorced from an ability to identify targets of 
%interest. Each new target is 'handed off' to the algorithm in its first
%detection frame. From that point on, it is the job of the algorithm to
%track it and to recognize when the object has left the frame. It receives
%no follow up information from the 'handoff' function other than new object
%declarations. Tracking is performed using a boosting framework. A variety
%of object properties are identified as linearly combined to form a
%stronger object association tool. Details follow in the code below.
%NOTE: Outer for loops in the program allow for implementation of more than
%one track algorithm in the same framework. A goal of this program was to
%have a high level of flexibility to incorporate new tracking approaches in
%the future. It will also allow for experimentation with combining tracking
%approaches.
%-------------------------------------------------------------------------%

%Bookkeeping:

%Add the working folder to the path 
%addpath('C:\Users\Pushpak\Documents\MATLAB\OTCBVS\Dataset 1');
%cd('C:\Users\Pushpak\Documents\MATLAB\OTCBVS\Dataset 1');

%Select the desired folder
c_folder = '00001';

%First access the Dataset 1/5 Directory list for the selected folder:
cd(c_folder);
d_list = dir('*.bmp');  %Directory listing of all frames
close all;
clearvars -except c_folder d_list;

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

%Create handoffs and ground truth matrices
GTPlots = GetGTPlots(c_folder);
Handoffs = Handoff(c_folder);

%Create the tracking lists: Multiple lists can be maintained. Each list
%would correspond to a specific tracking method
num_tmethods = 1;   %Define the number of tracking methods
t_list = zeros(0,0,num_tmethods);
background_status = 1;  %Sets background status (0 for bg, 1 for fg)
bgf_count = 1;          %Background frame count  
templates{1,1,3} = zeros(1);         %Stores image templates to extract properties
stats{1,1} = zeros(1);

disp('Program Start...');    
    
for i=1:numel(d_list)
    
    I_curr = double(imread(d_list(i).name));
    %---------------------------------------------------------------------%
    %Preprocessing: 
    %Performs background subtraction, binarization off of an empirically
    %defined threshold, erosion to remove small noise objects, and finally
    %displays the processed image as well as the original image in seperate
    %figures
    
    %As a default, assume objects are not detected (background frame)
    background_status = 1;
    
    %Take the difference between frame and bg to find foreground
    I_diff = abs(I_curr-I_bg);    
    %---------------------------------------------------------------------%
    
    
    %Step 1: For each new image, first determine if there are any new
    %handoffs. After handoff its the algorithms job to maintain track

    for n = 1:size(Handoffs,2)
        if Handoffs(i,n) == 1                
            n_size = size(t_list,1)+1;
            for x = 1:size(t_list,3)
                t_list(n_size,1,x) = GTPlots{i,n}.x(1);
                t_list(n_size,2,x) = GTPlots{i,n}.y(1);
                t_list(n_size,3,x) = GTPlots{i,n}.c(3);
                t_list(n_size,4,x) = GTPlots{i,n}.c(4);
                t_list(n_size,5,x) = GTPlots{i,n}.x(1);
                t_list(n_size,6,x) = GTPlots{i,n}.y(1);
                t_list(n_size,7,x) = GTPlots{i,n}.x(1);
                t_list(n_size,8,x) = GTPlots{i,n}.y(1);

                %All peripheral statistics must be determined for any new
                %detections at this point
                
                %First, create correlation templates:
                ulX = t_list(n_size,1,x);
                lrX = t_list(n_size,1,x) + t_list(n_size,3,x);
                ulY = t_list(n_size,2,x);
                lrY = t_list(n_size,2,x) + t_list(n_size,4,x);
                templates{n_size,x,1} = I_diff(ulY:lrY,ulX:lrX);
                templates{n_size,x,2} = I_curr(ulY:lrY,ulX:lrX);
                templates{n_size,x,3} = templates{n_size,x,1};
                
                %Accumulate stats identifying this object
                S = regionprops((templates{n_size,x,1}>45),'ConvexArea','Extent',...
               'MajorAxisLength','MinorAxisLength','Orientation','Perimeter');
                if isempty(S) == 0
                    if length(S) == 1
                        stats{n_size,x} = S;
                    else
                        %If there is more than one object found, select the one
                        %with the largest convex area and use its properties
                        maxi = S(1).ConvexArea;
                        ind = 1;
                        for v = 2:length(S)
                            if S(v).ConvexArea>maxi
                                ind = v;
                            end
                        end
                        stats{n_size,x} = S(ind);
                    end   
                end
            end
        end
    end
    figure(1),imshow(I_curr,[0 255]);
    hold on;
    %figure(2),imshow(I_diff,[-255 255]);
    figure(2),imshow(I_diff>45,[0 1]);
    hold on;    

    %---------------------------------------------------------------------%
    %Perform Tracking:
    %Track list updates for each tracking technique at this point
    
    thresh_l = [];
    min_l = zeros(size(t_list,1),1);
    %Identify new objects of interest!
    c_t = I_diff>45;
    c_c = regionprops(c_t, 'Centroid','ConvexArea','Extent',...
          'MajorAxisLength','MinorAxisLength','Orientation','Perimeter');
    val_thresh = .7;
    
    %Remove objects that are too small!
    for d = length(c_c):-1:1
        if isempty(c_c(d)) == 0   
              if c_c(d).ConvexArea < 10
                  c_c(d) = [];
              end
        end
    end
    
    for l = 1:size(t_list,3)
            for n=1:size(t_list,1)
            if t_list(n,1,l) ~= -1    
                %Find correlation peaks
                c = normxcorr2(templates{n,l,1}, I_diff);
                c_i = normxcorr2(templates{n,l,2}, I_curr);
                
                %Make the images of equal size:
                ys = round((-size(I_diff,1)+size(c,1))/2);
                xs = round((-size(I_diff,2)+size(c,2))/2);
                c = c(ys:ys+size(I_diff,1)-1,xs:xs+size(I_diff,2)-1);
                
                ys = round((-size(I_curr,1)+size(c_i,1))/2);
                xs = round((-size(I_curr,2)+size(c_i,2))/2);
                c_i = c_i(ys:ys+size(I_curr,1)-1,xs:xs+size(I_curr,2)-1);
                
%This was another method I considered using, but later disregarded. This
%method identified objects of interest based on correlation peaks but was
%less effective than later methods.
%                c_c = [];              
%                 for b = 1:length(c_cc)
%                     sx = t_list(n,3,l);
%                     sy = t_list(n,4,l);
%                     
%                     ys = max(1,round(c_cc(b).Centroid(2)-sy/2));
%                     xs = max(1,round(c_cc(b).Centroid(1)-sx/2));
%                     i1 = I_diff(ys:min(size(I_diff,1),round(ys+sy/2)),...
%                                 xs:min(size(I_diff,2),round(xs+sx/2)));
%                             
%                     c_n = regionprops(i1>45, 'Centroid','ConvexArea','Extent',...
%                     'MajorAxisLength','MinorAxisLength','Orientation','Perimeter');
%                    if isempty(c_n) == 0
%                        if length(c_n) == 1
%                            c_c{b} = c_n;
%                        else
%                            %If there is more than one object found, select the one
%                            %with the largest convex area and use its properties
%                            maxi = c_n(1).ConvexArea;
%                            ind = 1;
%                            for v = 2:length(c_n)
%                                if c_n(v).ConvexArea>maxi
%                                    ind = v;
%                                end
%                            end
%                            c_c{b} = c_n(ind);
%                        end
%                    end
%                    %figure(8), imshow(i1>45,[0 1]);
%                    %dd = waitforbuttonpress;
%                 end


                
                %For all detecctions peaks, also get other stats. Then use
                %the full range of stats to assign track correspondences
                [m0 ind1] = max(c);
                [m1 ind2] = max(m0);
                peakc = [ind1(ind2) ind2];
                
                [m2 ind3] = max(max(c_i));
                
                
                %Now determine for the object in question what if any
                %object in this frame is the best match:
                
                %C_C is the moderate correlation list. It is a list of starting
                %guesses
                thresh = [];                   
                x_predicted = t_list(n,7,l);
                y_predicted = t_list(n,8,l);                  
                aspect_predicted = stats{n,l}.MajorAxisLength/stats{n,l}.MinorAxisLength;
                
                %Boosting parameters: Set by guesses initially, and later
                %determined to minimize total error:
                %These parameters should be altered to desired values!
                a1 = .1; %a1 = .05;
                a2 = .1; %a2 = .05;
                a3 = 0; %a3 = .2;
                a4 = 0; %a4 = 0.0;
                a5 = 0.1; %a5 = .00;
                a6 = 0;  %a6 = 0.00;
                a7 = 0; %a7 = 0.00;
                
                %Assume that in the first four frames, background
                %subtracted data is inaccurate. After that point assume it
                %is more accurate. In this simplified form it has little
                %overall impact but could be expanded in the future to
                %accurately weight data based on its certainty. The instead
                %of discrete changes in the weights at the 5th frame, they
                %can change continuously over time.
                
                if i < 5
                    a8 = 0;
                    a11 = .9;
                    a9 = 0;
                    a10 = 0;
                else
                a8 = .45; %a8 = .1;                
                a11 = .45; %a11 = .55;                
                a9 = .1; %a9 = .025;
                a10 = .1; %a10 = .025;
                end

                
                for j = 1:length(c_c)
                if isempty(c_c(j)) == 0 
                    aspect_ratio = c_c(j).MajorAxisLength/c_c(j).MinorAxisLength;
                    
                    window = 4;
                    x_min = max(round((c_c(j).Centroid(1)))-window,1);
                    x_max = min(round((c_c(j).Centroid(1)))+window,size(c,2));
                    y_min = max(round(c_c(j).Centroid(2))-window,1);
                    y_max = min(round(c_c(j).Centroid(2))+window,size(c,1));
                    corr = max(max(c(y_min:y_max,x_min:x_max)));
                    corr_i = max(max(c_i(y_min:y_max,x_min:x_max)));
                    
                    %Initially started with 11 boosting parameters,
                    %calculate V!
                    thresh(j) = a1* abs(c_c(j).Centroid(1) - x_predicted)/size(I_diff,2) +...
                                a2 * abs(c_c(j).Centroid(2) - y_predicted)/size(I_diff,1) +...
                                a3 * abs(c_c(j).ConvexArea - stats{n,l}.ConvexArea)/max(.1,stats{n,l}.ConvexArea)+...
                                a4 * abs(c_c(j).Extent - stats{n,l}.Extent)/max(.1,stats{n,l}.Extent)+...
                                a5 * abs(aspect_ratio - aspect_predicted)/max(.1,aspect_predicted)+...
                                a6 * abs(c_c(j).Orientation - stats{n,l}.Orientation)/90+...
                                a7 * abs(c_c(j).Perimeter - stats{n,l}.Perimeter)/max(.1,stats{n,l}.Perimeter)+...
                                a8 * abs(1 - corr)/1+...
                                a9 * abs(c_c(j).Centroid(1) - peakc(1))/size(I_diff,2)+...
                                a10 * abs(c_c(j).Centroid(2) - peakc(2))/size(I_diff,1)+...
                                a11 * abs(1 - corr_i)/1;
                else
                    %This obect is now unselectable if we didn't see viable
                    %objects in the window:
                    thresh(j) = 10;
                end
                end
                if isempty(thresh) == 0
                    %Find the minimum association
                    [val t_ind] = min(thresh);
                    thresh_l{n} = thresh;
                    min_l(n) = val;
                    
                    t_list(n,1,l) = min(size(I_diff,2),max(1,round(c_c(t_ind).Centroid(1)-t_list(n,3,l)/2)));
                    t_list(n,2,l) = min(size(I_diff,1),max(1,round(c_c(t_ind).Centroid(2)-t_list(n,4,l)/2)));
                    if val > val_thresh
                        %If this threshold is reached, we could not find an
                        %associate track point. It is then assumed that this
                        %point has now left the screen. If this is the case, we
                        %delete it by making the row full of -1's. This
                        %indicates to the  program to ignore it in the future:
                        t_list(n,:,l) = -1 * ones(1,8);
                        min_l(n) = 10;
                    end
                else
                        %If this threshold is reached, we could not find an
                        %associate track point. It is then assumed that this
                        %point has now left the screen. If this is the case, we
                        %delete it by making the row full of -1's. This
                        %indicates to the  program to ignore it in the
                        %future:
                        t_list(n,:,l) = -1 * ones(1,8);
                        %Fill min list with an arbitrary large number
                        min_l(n) = 10;
                end
                
            end
            end
            
            %Conflict = 0 indicates that there are no repeat associations.
            %Conflict = 1 indicates that a repeat association needs to be
            %reconciled
            conflict = 1;
            while conflict == 1
                conflict = 0;
                for n = 1:length(thresh_l)
                    if isempty(thresh_l{n}) == 0 && min_l(n) <10
                        for k=1:length(thresh_l)
                            if isempty(thresh_l(k)) == 0 && min_l(k) <10
                            if n ~= k
                            %Stat compilation now makes it possible to evaluate which
                            %is the most likely associated object

                            %Centroids must not be within 10 pixels or we have
                            %an assumed conflict
                            t_dis = 10;
                            if abs(t_list(n,1,l)-t_list(k,1,l))<=t_dis && ...
                               abs(t_list(n,2,l)-t_list(k,2,l))<=t_dis
                                
                                %Then there is a conflict! Sort out the
                                %conflict by giving the association to the
                                %track_object that has the lower val with
                                %that image patch
                                if min_l(n) > min_l(k)
                                    [val t_ind] = min(thresh_l{n});
                                    thresh_l{n}(t_ind) = 10;
                                    [val t_ind] = min(thresh_l{n});                   
                                    t_list(n,1,l) = min(size(I_diff,2),...
                                        max(1,round(c_c(t_ind).Centroid(1)-t_list(n,3,l)/2)));
                                    t_list(n,2,l) = min(size(I_diff,1),...
                                        max(1,round(c_c(t_ind).Centroid(2)-t_list(n,4,l)/2)));
                                    %Update the new min and make sure it
                                    %still has a viable association
                                    if val > val_thresh
                                        t_list(n,:,l) = -1 * ones(1,8);
                                        min_l(n) = 10;
                                    else
                                        min_l(n) = val;
                                    end
                                    conflict = 1;
                                end
                            end
                            end
                            end
                        end
                    else
                    min_l(n) = 10;
                    end
                end
            end
    end         
disp('NEXT');
    %---------------------------------------------------------------------%
    %Display all tracks: At this point tracks have been maintained 
    %seperately in the t_list matrix. Also ground truths have been
    %maintained. This section displays all these tracks simultaneously on
    %the frame of interest.
    p = 1;
    while p <= size(GTPlots,2)
        if isempty(GTPlots{i,p}) == 0
            figure(1), plot(GTPlots{i,p}.x,GTPlots{i,p}.y, 'r');
            hold on;
        end
        p = p+1;
    end
        
    for x = 1:size(t_list,3)
    %Plot the elements of each t_list here
        p = 1;
        while p <= size(t_list,1)
            %Do NOT plot if the row is empty or indicated as a deleted
            %track
            if sum(t_list(p,:,x)) ~= 0 && t_list(p,1,x) ~= -1
                ulX = t_list(p,1,x);
                lrX = t_list(p,1,x) + t_list(p,3,x);
                ulY = t_list(p,2,x);
                lrY = t_list(p,2,x) + t_list(p,4,x);
                figure(1), plot([ulX lrX lrX ulX ulX],[ulY ulY lrY lrY ulY],'Color', [.5,.1,1/x]);
                hold on;
            end
            p = p+1;
        end
        %At this point update position prediction for the next frame
        t_list(:,:,x) = tlupdate(t_list(:,:,x),I_diff);
    end
    
    
    %---------------------------------------------------------------------%
    %Update All Track Statistics based on this frame:
    %First update the correlation templates
    
    for l = 1:size(t_list,3)
        for n = 1:size(t_list,1)
        if t_list(n,1,l) ~= -1 
            %First create image templates for each object
            ulX = max(1,t_list(n,1,l));
            lrX = min(size(I_diff,2),t_list(n,1,l) + t_list(n,3,l));
            ulY = max(t_list(n,2,l),1);
            lrY = min(size(I_diff,1),t_list(n,2,l) + t_list(n,4,l));
            templates{n,l,1} = I_diff(ulY:lrY,ulX:lrX);
            templates{n,l,2} = I_curr(ulY:lrY,ulX:lrX);

            if isempty(templates{n,l,3}) == 0
                templates{n,l,3} = (templates{n,l,1} + templates{n,l,2})/2;
            else
                templates{n,l,3} = templates{n,l,1};
            end
            
            %Accumulate stats identifying this object
            S = regionprops((templates{n,l,1}>45),'ConvexArea','Extent',...
           'MajorAxisLength','MinorAxisLength','Orientation','Perimeter');
            
           if isempty(S) == 0
                if length(S) == 1
                    stats{n,l} = S;
                else
                    %If there is more than one object found, select the one
                    %with the largest convex area and use its properties
                    maxi = S(1).ConvexArea;
                    ind = 1;
                    for v = 2:length(S)
                        if S(v).ConvexArea>maxi
                            ind = v;
                        end
                    end
                    stats{n,l} = S(ind);
                end  
           end
        end
        end
    end
    %Now all objects in the track list have 6 identified properties, as
    %well as an image template, an expected position in the next frame.
    %These characteristics are used for tracking.

    %---------------------------------------------------------------------%
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
    %---------------------------------------------------------------------%

    %Pause to allow the user to view the frame. Proceed on a button press:
    k = waitforbuttonpress; 

end


            
