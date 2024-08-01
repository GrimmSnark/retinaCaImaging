function [riseTime, decayTime, riseThreshIndxRelative, fallThreshIndx] = riseTimeDecayFinderTriangle(riseSample, fallSample, rate)
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
%          riseThreshIndxRelative - rise start index relative to the peak
%                                   index
%
%          fallThreshIndx - decay stop index relative to the peak index

%% data processing

% smooth data
riseSampleSm = smoothdata(riseSample);
fallSampleSm = smoothdata(fallSample);

% get the elbow points
riseThreshIndx = triangle_threshold(riseSampleSm, 'L',0)-1;
% [~, fallThreshIndx] = min(fallSampleSm);
fallThreshIndx = triangle_threshold(fallSampleSm, 'R',0)+5;


% make riseThresholdINdx relative to the peak index
riseThreshIndxRelative = (length(riseSample)- riseThreshIndx);


%% uncomment for plotting
% subplot(211)
% plot(riseSampleSm)
% hold on
% scatter(riseThreshIndx,riseSampleSm(riseThreshIndx));
% hold off
% 
% subplot(212)
% plot(fallSampleSm)
% hold on
% plot(gradient(fallSampleSm))
% scatter(fallThreshIndx,fallSampleSm(fallThreshIndx));
% hold off

% work out actual times etc
riseTime = riseThreshIndxRelative/ rate; % makes the rise time relative to the peak point NOT the total index length!!!
decayTime = fallThreshIndx/rate;

end