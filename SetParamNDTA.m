%SetParmNDTA: Optimization algorithm for NoDetect_Track
%Methodology: 
%SetParmNDTA iterates rapidly through a wide range of input parameters
%(a1-a11) which represent weights of boosting parameters, val_thresh, and
%similar parameters in NoDetect_Track. The automated version of 
%NoDetect_Track is called for performance estimation. It is run on only
%image sequence 1 due to the already extensive computation time required.
%Error in number of objects detected in each frame is calculated for
%each scenario (and parameter input combination). Ground truth data is then
%extracted and compared to the calculated data. The minimum error input 
%parameter combinations are estimated and shown below. This function serves 
%as a simplified brute force parameter optimization approach.

%NOTE: The full nested loop algorithm below was not run due to time
%constraints. Instead a wide variety of combinations were run and the
%resulting data were analyzed to get a sense of the optimal parameters. In
%general 4 or 5 loops would be run, and optimal parameters for those cases
%would be estiamted. Then I would iterate back and run differnt loop
%combinations. The resulting estimated optimization is a first pass at
%finding optimal weighting values. A more efficient optimization
%methodology will be implemented in the future. Commented in and out loops
%are left as they were on the last optimization run I peformed. In general
%for each useful optimization run, I would store the datapoints (generally
%approximately 300 - 500) in a matrix V_NUMBER. Then in later effective
%runs I would either write over unneeded data or increment to V2 and V3 and
%so on...

%Add the working folder to the path --> Change to appropriate
addpath('C:\Users\Pushpak\Documents\MATLAB\OTCBVS\Dataset 1');

V8 = [];
%Extract ground truth data for image sequence 1:        %26 could be 3 or 4
N_T = [5 5 1 1 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 4 3 4 7 7 7 7];
c_folder = '00001';
count = 1;

%for a1 = 0:.1:.1
    a1 = .1;
    a2 = a1;
        for a3 = 0:.1:.3
        %a3 = 0;
            %for a4 = 0:.1:.1
            a4 = 0;
                %for a5 = 0:.1:.1
                a5 = .05;
                    %for a6 = 0:.05:.1
                    a6 = 0;
                        %for a7 = 0:.0:.1
                            a7 = 0;
                            for a8 = 0.4:.05:.5
                            a8 = .45;
                                for a9 = 0:.1:.1
                                    %a9 = .1;
                                    a10 = a9;
                                    for a11 = .4:.05:.55
                                    %a11 = .45;
                                        for val_thresh = .5:.1:.9
                                            tdiff = 45;
                                            
                                            GTPlots = GetGTPlots(c_folder);
                                            Handoffs = Handoff(c_folder);
                                            O = NoDetect_Track_A(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,val_thresh,GTPlots,Handoffs, tdiff);
                                            %4 or 3 is acceptable for Frame
                                            %26 of the video:
                                            if O(26) == 4
                                                O(26) = 3;
                                            end
                                            V8(count,:) = [sum(abs(N_T - O)) a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 val_thresh tdiff];
                                            disp(count)
                                            count = count + 1;
                                            
                                        end
                                    end
                                end
                            end
                        %end
                    %end
                %end
            %end
        end
%end

%Optimally Determined Values:
%[0.1,0.1,0,0,0.05,0,0,0.45,0.1,0.1,0.45,0.7,45;]
                                  
