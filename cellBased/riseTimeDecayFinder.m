function [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinder(riseSample, fallSample, rate)
% Function tries to find rise time and tau decay for calcium events.
% Heavily based on the risetime and falltime function. 
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


%% lower threshold rise value crossing

% get lower/upper state values
[~,~,~,~,~, stateLevelsRise] = risetimeState(riseSample);

%%%% take from plotFinalize
amp = diff(stateLevelsRise);
thresholdLevel = stateLevelsRise(1) + 2 * amp / 100;
%%%%

% find index for when the state level is crossed
crossIndxRise = find(riseSample>thresholdLevel,1);

riseTimeIndx = length(riseSample)- crossIndxRise;

riseTime = riseTimeIndx / rate; % rise time in seconds

%% lower threshold fall value crossing

% get lower/upper state values
[~,~,~,~,~, stateLevelsFall] = falltimeState(fallSample);

%%%% take from plotFinalize
amp = diff(stateLevelsFall);
thresholdLevel = stateLevelsFall(1) + 2 * amp / 100;
%%%%

% find index for when the state level is crossed
crossIndxFall = find(fallSample<thresholdLevel,1);

decayTime = crossIndxFall / rate; % decay time in seconds

end