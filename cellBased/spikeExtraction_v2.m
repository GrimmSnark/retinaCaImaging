function spikeExtraction_v2(exStructPath)

% This function completes 'spike' detection on calcium activity of
% idenitified cells. Calculates spike amplitudes, full width at half max,
% tau rise and decays
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs-  exStructPath: filepath for exStruct.mat to process

%% defaults

zscoreLim = 12;
% sIQRLim = 0.04;
plotOffset = 0.2;


%% load exStruct

if nargin < 1 || isempty(exStructPath)
    [file, path] = uigetfile({'*.mat'},...
        'Image File Selector');

    exStructPath = fullfile(path,file);
end

exStruct = load(exStructPath);
exStruct = exStruct.exStruct;
exStruct.cells.zscoreLim = zscoreLim;

%% analyse each cell
exStruct.spikes = [];

for c= 1:exStruct.cellCount
    dFTrace = exStruct.cells.dF(c,:);
    dFTrace = sgolayfilt(dFTrace,5,exStruct.rate+1);

    %% zscore
    [zScore(c,1), sIQR(c,1) ]= dF_zscore(dFTrace);
    exStruct.cells.zScore = zScore;
    exStruct.cells.zScoreThresholded = zScore > zscoreLim;


    %% 'spike' detection

    dF_SD = std(dFTrace);
    % spike detection
    [spikeAmp{c},spikeLocs{c},spikeWidths{c}] = findpeaks(dFTrace, "MinPeakProminence",dF_SD/2);


    % get tau rise time and decay
    % if ~isempty(spikeLocs{c}) && zScore(c) > zscoreLim
    if ~isempty(spikeLocs{c}) && zScore(c)
        
        % for each spike
        for sp = 1:length(spikeLocs{c})

            sampleSt = spikeLocs{c}(sp)- (exStruct.rate * 1); % take 1 sec before peak onset
            spikeIndex = exStruct.rate * 1;

            % corection if peak is less than a second before start
            if sampleSt <= 0
                sampleSt = 1;
            end

            sampleEnd = spikeLocs{c}(sp)+ (exStruct.rate * 3); % take 4 sec after peak onset

            % corection if sample end is greater than dFtrace end
            if sampleEnd > length(dFTrace)
                sampleEnd = length(dFTrace);
            end

            % get single transient trace
            spikeTrace = dFTrace(sampleSt:sampleEnd);

            % rise/decay time threshold
            [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinderTriangle(spikeTrace, spikeIndex, exStruct.rate);


            % add into structure
            riseTimeStruct{c}(sp) = riseTime;
            decayTimeStruct{c}(sp) = decayTime;
            riseTimeIndxStruct{c}(sp) = spikeLocs{c}(sp)-riseTimeIndx;

            % check if last decay is the end of trace
            if (spikeLocs{c}(sp)+crossIndxFall) > length(dFTrace)
                decayTimeIndxStruct{c}(sp) = length(dFTrace);
            else
                decayTimeIndxStruct{c}(sp) = spikeLocs{c}(sp)+crossIndxFall;
            end
        end
    else
        riseTimeStruct{c} = [];
        decayTimeStruct{c} = [];
        riseTimeIndxStruct{c} = [];
        decayTimeIndxStruct{c} = [];
    end

    % add firing rate

    if any(riseTimeStruct{c})
        firingRate(c) = length(riseTimeStruct{c})/ (length(exStruct.xyShifts)/exStruct.rate);
        meanAmp(c) = sum(spikeAmp{c})/length(spikeAmp{c});
    else
        firingRate(c) = 0;
        meanAmp(c) = 0;
    end
end

%%  clean decay structure so that it does not exceed the next rise time....
for i = 1:length(decayTimeIndxStruct)
    for x = 1:length(decayTimeIndxStruct{i})-1
        if decayTimeIndxStruct{i}(x) > riseTimeIndxStruct{i}(x+1)
            decayTimeIndxStruct{i}(x) = riseTimeIndxStruct{i}(x+1)-1;
        end
    end
end

% recalucate the decay time structure
res=cellfun(@(x,y) y-x,spikeLocs,decayTimeIndxStruct, 'UniformOutput',false);
decayTimeStruct = cellfun(@(x) x/exStruct.rate,res, 'UniformOutput',false);

%% add everything into exStruct
exStruct.spikes = [];

if isfield(exStruct,  'spikesSmall')
    exStruct = rmfield(exStruct, 'spikesSmall');
end

% large spikes
exStruct.spikes.spikeAmp = spikeAmp;
exStruct.spikes.spikeLocs = spikeLocs;
exStruct.spikes.spikeWidths = spikeWidths;
exStruct.spikes.riseTime = riseTimeStruct;
exStruct.spikes.decayTime = decayTimeStruct;
exStruct.spikes.riseTimeIndx = riseTimeIndxStruct;
exStruct.spikes.decayTimeIndx = decayTimeIndxStruct;
exStruct.spikes.firingRate = firingRate; % firing rate in spk/s
exStruct.spikes.meanSpikeAmp = meanAmp; % mean spike dF amplitude

%% Do the coherence metrics
exStruct = waveCellCoherence(exStruct);


%% plot everything
zscoreFig = figure('units','normalized','outerposition',[0 0 1 1]);
zscoreAx = gca;
hold on
axis tight
%% plot based on zcore
for c= 1:exStruct.cellCount

    dFTracePlot = exStruct.cells.dF(c,:) + (c-1)*plotOffset;
    offsetVals(c) = (c-1)*plotOffset;

    if zScore(c) > zscoreLim
        plot(zscoreAx,dFTracePlot)

        text(zscoreAx,50,dFTracePlot(50),0,num2str(zScore(c)));

        try
            % large spikes
            plot(zscoreAx,spikeLocs{c}(:), dFTracePlot(spikeLocs{c}(:)), 'r*');

            plot(zscoreAx,riseTimeIndxStruct{c} , dFTracePlot(riseTimeIndxStruct{c}), 'b^');
            plot(zscoreAx,decayTimeIndxStruct{c} , dFTracePlot(decayTimeIndxStruct{c}), 'gv');

        catch
        end
    else
        plot(zscoreAx,dFTracePlot, '--','Color', [0.5 0.5 0.5]);
        text(zscoreAx,50,dFTracePlot(50),0,num2str(zScore(c)));
    end
end

%% plot waves

if sum(exStruct.cells.zScoreThresholded) > 0
    spikesSorted = exStruct.wavesMetrics.spikesSorted;

    for d = 1:max(spikesSorted(:,3))
        wavePoints = spikesSorted(spikesSorted(:,3)==d,:);

        cellIndex = wavePoints(:,2);% to add amplitude offset to align to the top of the peak

        wavePoints(:,2) =  ((wavePoints(:,2))*plotOffset)-plotOffset;
        wavePoints(:,3) = [];

        % adding peak Y offset to better align
        for xx = 1: length(cellIndex)
            wavePoints(xx,2) = wavePoints(xx,2) + spikeAmp{cellIndex(xx)}(find(spikeLocs{cellIndex(xx)}==wavePoints(xx,1)));
        end

        wavePoints = sortrows(wavePoints,2);

        plot(wavePoints(:,1)',wavePoints(:,2)', 'LineStyle',':');
    end
end

tightfig;

fileStruct = dir(exStruct.filePath);
saveas(zscoreFig,fullfile(fileStruct.folder,[fileStruct.name(1:end-4) 'DF_fig.tif']));

%% save data
save(exStructPath, "exStruct", '-v7.3');

end