%Track.m:
%This was a first pass at meanshift tracking and also correlation tracking.
%However although the correlation (shown, meanshift was subsequently
%removed and called in single_track for better demonstration) was able to
%identify points at least somewhat, there was no logical maintenance of a
%track list included. This development work spurred the creation of single
%track to include uncertainty and a methodology for adding, removing, and
%association objects on the track list.

    SEQ_DIR = '00001';

    clear T;
    clear S;
    
    %Grab the Groundtruth file with associated sequence
    filename = fullfile(SEQ_DIR, 'groundTruth.txt');
    fid = fopen(filename);
    
    %Cycle through comments in the file:
    line = fgets(fid);
    while line(1) == '%'
    line = fgets(fid);
    end
    
    %Read the number of images in the video
    numImages = sscanf(line, '%d', 1);
    start = 5;
    
    for i =start:10;
    %for i=1:numImages
        
        %Load the image name and number of boxes
        imageName = fscanf(fid, '%c',13);
        numBoxes = fscanf(fid, '%d', 1);
        
        %Display the image
        fname = fullfile(SEQ_DIR, imageName);
        Im = imread(fname);
        imagesc(Im);
        colormap('gray');
        axis('off');
        title(sprintf('Image %d', i));
        hold on;
        
        %Load the ground truth boxes
        for j=1:numBoxes
            tmp = fscanf(fid, '%c',2); %% [space](
            coords = fscanf(fid, '%d %d %d %d');
            tmp = fscanf(fid, '%c',1); %% )
            ulX=coords(1); ulY=coords(2);
            lrX=coords(3); lrY=coords(4);
            boxes{j}.X = [ulX lrX lrX ulX ulX]';
            boxes{j}.Y = [ulY ulY lrY lrY ulY]';
        end
        
        tmp = fgetl(fid); %% get until end of line
        
        %Display the ground truth boxes (Colored RED)
        for j=1:numBoxes
            plot(boxes{j}.X, boxes{j}.Y, 'r');
            
            %Create Tracks - for multiple object tracking
            
            %If there are no tracks to begin with, then all detected
            %objects should be added as tracks
            %if(size(T,3)>0)
            if(i == start)
                %Extract template based on groud truth data
                Tem = Im(boxes{j}.Y(1):boxes{j}.Y(3),...
                         boxes{j}.X(1):boxes{j}.X(3));
                T{j} = Tem;                
                
                S(j,2) = boxes{j}.X(3) - boxes{j}.X(1);
                S(j,1) = boxes{j}.Y(3) - boxes{j}.Y(1);
            end
        end
        
        %Apply correlation tracking - disregard tracks with low correlation
        %--> The track object has probably left the screen
        t_s = 0; 
        
        if(i>start && isempty(T)== 0)
            for n=1:size(T,2)    
                g = normxcorr2(T{1,n - t_s}, double(Im));
                
                if (max(max(g))>.70)
                    gT = g == max(g(:));
                    [iy ix] = find(gT == 1);
                    numel(ix);
                
                    sx = S(n-t_s,2);
                    sy = S(n-t_s,1);
                    
                    plot([ix ix ix-sx ix-sx ix],...
                         [iy iy-sy iy-sy iy iy], 'y');                   
                    %Update the template based on current results
                    T{1,n-t_s} = (T{1,n-t_s} + Im(iy:iy+sy,ix:ix+sx))/2;
                else
                    T(:,n - t_s) = [];
                    S(n - t_s,:) = [];
                    t_s = t_s + 1;
                end
            end
        end
         
        %Pause to set frame rate
        pause(.5);
        hold off;
        clear Im;
        
    end

    fclose(fid);
%end
