%Meanshift Function:
%Created a simple convergence algorithm which looks in windows of varying
%sizes. This function is called by wrapper programs such as single
%track.This algorithm identifies the first window that is stable.

function [s_win x_c y_c] = camshift(I_curr, s_win)

%Constant Parameters:
T = .5;               % Threshold of convergence (in pixels)
del = T;             %Pixel delta, defaults at T
max_iterations = 200;

%Input the initial search window (location and size)
ix = s_win(1);      %Initial x location
iy = s_win(2);      %Initial y location
sx = s_win(3);      %Search window size in x
sy = s_win(4);      %Search window size in y

for x_f = .5:.1:1.5
for y_f = .5:.1:1.5

    x_c = round(ix + sx/2); 
    y_c = round(iy + sy/2);

    iterations = 0;

    %While the solution is still converging, as long as we haven't been
    %searching for too many iterations
    while (iterations < max_iterations && del >= T)
        
        %Save the previous ix and iy:
		x_p = x_c;
        y_p = y_c;
        
        % Compute centroid of search window
        x_min = round(x_c - x_f * sx); 
        x_max = round(x_c + x_f * sx);
        y_min = round(y_c - y_f * sy);
        y_max = round(y_c + y_f * sy);
        
		TS = double(0);
		for i = x_min:x_max
            for j = y_min:y_max
                if i < size(I_curr,2)&& j < size(I_curr,1)&& i > 1&& j >  1
                    %Calculate the sum over the area:
                    TS = TS + double(I_curr(j,i));
                end
            end
		end
		
		I_x = double(0);
		for i = x_min:x_max
            for j = y_min:y_max
                if i < size(I_curr,2)&& j < size(I_curr,1)&& i > 1&& j >  1
                    %Calculate the weighted value:
                    I_x = I_x + i * double(I_curr(j,i));
                end
            end
		end
		
        I_y = double(0);
        for i = x_min:x_max
            for j = y_min:y_max
                if i < size(I_curr,2)&& j < size(I_curr,1)&& i > 1&& j >  1
                    %Calculate the weighted value:
                    I_y = I_y + j * double(I_curr(j,i));
                end
            end
        end
		
        %Find the centroid:
		x_c = round(I_x/TS);
		y_c = round(I_y/TS);
        
       
        %Calculate Current Error
        del = abs(x_p-x_c) + abs(y_p-y_c);
        iterations = iterations + 1;
    end
    
    if del <= T
        break;
    end
        
end
    if del <= T
        break;
    end
end

sx = min(max(round(x_f * sx),1),size(I_curr,2));
sy = min(max(round(y_f * sy),1),size(I_curr,1));
ix = min(max(round(x_c - x_f * sx/2),1),size(I_curr,2));
iy = min(max(round(y_c - y_f * sy/2),1),size(I_curr,1));

%Set return varaibles for new window location
s_win(1) = ix;      %Initial x location
s_win(2) = iy;      %Initial y location
s_win(3) = sx;      %Search window size in x
s_win(4) = sy;      %Search window size in y
end
    
    
       
	
    
