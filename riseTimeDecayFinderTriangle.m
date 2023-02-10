function [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinderTriangle(riseSample, fallSample, rate)
% Function tries to find rise time and tau decay for calcium events.
% Heavily based on triangle thresholding. 
% https://www.mathworks.com/matlabcentral/answers/250257-find-turning-point-in-data
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
% 
% Inputs: riseSample - calcium trace leading up to spike peak
%
%         fallSample - calcium trace from peak to decay afterwards
%
%         rate - imaging rate in Hz
%
% Outputs: riseTime - spike rise time in seconds
%
%          decayTime - spike decay time in seconds
%
%          riseTimeIndx - rise start index
%
%          crossIndxFall - decay stop index

%% data processing

% smooth data
riseSampleSm = smoothdata(riseSample);
fallSampleSm = smoothdata(fallSample);

% get the elbow points
riseThreshIndx = triangle_threshold(riseSampleSm, 'L',0);
fallThreshIndx = triangle_threshold(fallSample, 'R',0);


% work out actual times etc
riseTime = riseThreshIndx/ rate;
decayTime = fallThreshIndx/rate;

end