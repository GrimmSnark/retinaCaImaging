function [riseTime, decayTime, riseThreshIndxRelative, fallThreshIndx] = riseTimeDecayFinderTriangle(sample, peakIndex, rate)
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
% riseSampleSm = smoothdata(sample);
% fallSampleSm = smoothdata(fallSample);

% get the elbow points
riseThreshIndx = triangle_threshold(sample(1:peakIndex), 'L',0)-1;

if riseThreshIndx ==0
    riseThreshIndx = 1;
end

% smooth fall sample
fallSample = smoothdata(sample(peakIndex:end),1,"gaussian",3);
ft = fittype('a*exp(-b*t) + c','indep','t');

warning("off", 'all');
% try some exp fits
try
    [f, g] = fit( [1:length(fallSample)]', fallSample', ft);
catch
    [f, g] = fit( [1:length(fallSample)]', fallSample', 'exp1');
end

% if bad fit cut the fallsample in two (bad fit is ususally two peaks in
% the sample)
if g.rsquare < 0.7
    newInd = round(length(fallSample)/2);
    [f, g] = fit( [1:newInd]', fallSample(1:newInd)', 'exp1');
end

warning("on", 'all');

% no idea why this works....but it does
if f.b < 0
tt = round(-1/f.b) +5;
else
tt = round(1/f.b) +5;
end

% subplot(211)
% plot(f,[1:length(fallSample)], fallSample)
% title(gca, [ ' R-Squared: ' num2str(g.rsquare)] )

if isempty(tt)
    disp('Could not find local minima, defauting to end of peak trace')
tt = triangle_threshold(sample(peakIndex:end), 'R',0)+5;
end

fallThreshIndx = tt+peakIndex;

% make riseThresholdINdx relative to the peak index
riseThreshIndxRelative = (length(1:peakIndex)- riseThreshIndx);


%% uncomment for plotting 
% subplot(212);
% plot(sample)
% hold on
% scatter(riseThreshIndx,sample(riseThreshIndx));
% try
%     scatter(tt+peakIndex,sample(tt+peakIndex));
% catch
%     scatter(length(sample),sample(end));
% end
% hold off

% work out actual times etc
riseTime = riseThreshIndxRelative/ rate; % makes the rise time relative to the peak point NOT the total index length!!!
decayTime = fallThreshIndx/rate;

end