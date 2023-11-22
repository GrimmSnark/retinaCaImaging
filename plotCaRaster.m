function  [xPoints, yPoints] = plotCaRaster(spikes, col)


vertSpikeHeight = 0.5;
vertSpikePosition = 0;
%% Binary spike train matrix case. Initialize variables and set axes.
nTrials = size(spikes,1);
nTimes = size(spikes,2);

% Note: xlim and ylim are much, much faster than axis or set(gca,...).
xlim([0 nTimes+1]);
ylim([0 nTrials+1]);

%% Vertical Lines
% Find the trial (yPoints) and timebin (xPoints) of each spike
[trials,timebins] = find(spikes);
trials = trials';
timebins = timebins';
halfSpikeHeight = vertSpikeHeight/2;

xPoints = [ timebins;
    timebins;
    NaN(size(timebins)) ];
yPoints = [ trials - halfSpikeHeight + vertSpikePosition;
    trials + halfSpikeHeight + vertSpikePosition;
    NaN(size(trials)) ];

xPoints = xPoints(:);
yPoints = yPoints(:);
plot(xPoints,yPoints, 'Color', col);
set(gca,'YDir','reverse');

end