%tlupdate.m
%Methodology: Written as a general track list linear motion predictor
%function. This function is called by bother NoDetect_Track and
%Single_Track to perform estimation of object locations in subsequent
%frames. If necessary it could be expanded in the future to do non-linear
%motion prediction based on an established motion profile.

function [t_list] = tlupdate(t_list,I_diff)
%Updates the tracklist using linear prediction:
%Finds the delta between the previous frames, and adds this to the previous
%frame value in each dimension to predict its next location

    for i = 1:size(t_list,1)
        t_list(i,7) = max(1,min(size(I_diff,2),(t_list(i,1) - t_list(i,5))/1 + t_list(i,1)));
        t_list(i,8) = max(1,min(size(I_diff,1),(t_list(i,2) - t_list(i,6))/1 + t_list(i,2)));
    end

end

