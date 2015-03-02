%SetParm: Optimization algorithm for Single_Track
%Methodology: 
%SetParm iterates rapidly through a wide range of input parameters (p1 -
%p6) which relate to sizes of morphological processing windows and other
%similar parameters in Single_Track. The Single_Track called is an
%automated version. It is run on all 10 image sequences, and the frame
%error and distance error are calculated for each scenario (and parameter
%input combination). Ground truth data is then extracted and compared to
%the calculated data. The minimum error input parameter combinations are
%estimated and shown below. This function serves as a simplified brute
%force parameter optimization approach.

%RUN SINGLE TRACK OVER A WIDE PARAMETER RANGE
count = 1;
ALL = [];
p = [];

for p1 = 85:-10:35
    for p2 = 1:1:3
        for p3 = 1:1:3
            for p4 = 2:2:14
                for p5 = 2:2:8
                    for p6 = 0:1
                        

                            [ALL{count,1} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00001');
                            [ALL{count,2} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00002');
                            [ALL{count,3} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00003');
                            [ALL{count,4} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00004');
                            [ALL{count,5} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00005');
                            [ALL{count,6} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00006');                           
                            [ALL{count,7} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00007');
                            [ALL{count,8} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00008');                   
                            [ALL{count,9} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00009');
                            [ALL{count,10} p{count}] = Single_Track(p1,p2,p3,p4,p5,p6,'00010');
                            count = count + 1;                               
                            
                    end
                end
            end
        end
    end
end
    
%EXTRACT REQUIRED GROUND TRUTH DATA:
    GTPlots = [];

    GTPlots{1,1} = GetGTPlots1('00001');
    GTPlots{1,1}(1,:) = [];
    GTPlots{1,2} = GetGTPlots1('00002');
    GTPlots{1,2}(1,:) = [];
    GTPlots{1,3} = GetGTPlots1('00003');
    GTPlots{1,3}(1,:) = [];
    GTPlots{1,4} = GetGTPlots1('00004');
    GTPlots{1,4}(1,:) = [];
    GTPlots{1,5} = GetGTPlots1('00005');
    GTPlots{1,5}(1,:) = [];
    GTPlots{1,6} = GetGTPlots1('00006');
    GTPlots{1,6}(1,:) = [];
    GTPlots{1,7} = GetGTPlots1('00007');
    GTPlots{1,7}(1,:) = [];
    GTPlots{1,8} = GetGTPlots1('00008');
    GTPlots{1,8}(1,:) = [];
    GTPlots{1,9} = GetGTPlots1('00009');
    GTPlots{1,9}(1,:) = [];
    GTPlots{1,10} = GetGTPlots1('00010');
    GTPlots{1,10}(1,:) = [];
    
%COMPARE GROUND TRUTH DATA TO CALCULATED DATA:
TE = zeros(size(ALL,1),10);
M_Num = zeros(size(ALL,1),10);

for z = 1:size(ALL,1)
    for n = 1:10
        for x = 1:size(GTPlots{1,n}, 1)
            for y = 1:size(GTPlots{1,n}, 2)
                if isempty(GTPlots{1,n}{x,y}) == 0
                    if isempty(ALL{z,n}{1,x}) == 0
                        Err = abs(GTPlots{1,n}{x,y}.c(1) - ALL{z,n}{1,x}(:,1))+ ...
                              abs(GTPlots{1,n}{x,y}.c(2) - ALL{z,n}{1,x}(:,1));
                        TE(z,n) = min(Err) + TE(z,n);
                        m = length(ALL{z,n}{1,x}(:,1));
                    else
                        m = 0;
                    end
                else
                    break;
                end   
            end
            M_Num(z,n) = M_Num(z,n)+abs(y-m);
        end
    end
end

for g = 1:size(TE,1)
    rows(g) = sum(TE(g,:));
end

[min1 idx] = min(rows);
disp(p{idx});


for g = 1:size(M_Num,1)
    rows(g) = sum(M_Num(g,:));
end

[min2 idx] = min(rows);
disp(p{idx});


% The parameter output was:  [85     3     3     2     2     1]
%                            [65     1     1    12     8     0]

