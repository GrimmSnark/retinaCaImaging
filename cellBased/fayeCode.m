function fayeCode(exStructPath)
% This function completes 'spike' detection on calcium activity of
% idenitified cells. Calculates spike amplitudes, full width at half max,
% tau rise and decays
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs-  exStructPath: filepath for exStruct.mat to process

%% defaults

zscoreLim = 20;
% sIQRLim = 0.04;
frameDiffLimit = 10; % frame closeness limit for small and big spikes

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

    %% zscore
    [zScore(c,1), sIQR(c,1) ]= dF_zscore(dFTrace);

    %% 'spike' detection
    % min peak height
    %     dFTrace = movmedian(dFTrace,3);

    dF_SD = std(dFTrace);
    dFMean = mean(dFTrace);
    minHeight = dFMean + (dF_SD *3);
    minProm = 2 * dF_SD;

    % spike detection
    %     [spikeAmp{c},spikeLocs{c},spikeWidths{c}] = findpeaks(dFTrace, "MinPeakProminence",minProm);
    [spikeAmp{c},spikeLocs{c},spikeWidths{c}] = findpeaks(dFTrace,"MinPeakHeight",0.1, "MinPeakProminence",minProm);

    % get tau rise time and decay
    if ~isempty(spikeLocs{c})
        % for each spike
        for sp = 1:length(spikeLocs{c})

            sampleSt = spikeLocs{c}(sp)- (exStruct.rate * 1); % take 1 sec before peak onset

            % corection if peak is less than a second before start
            if sampleSt <= 0
                sampleSt = 1;
            end

            sampleEnd = spikeLocs{c}(sp)+ (exStruct.rate * 4); % take 4 sec after peak onset

            % corection if sample end is greater than dFtrace end
            if sampleEnd > length(dFTrace)
                sampleEnd = length(dFTrace);
            end

            % get rise time
            riseSample = dFTrace(sampleSt:spikeLocs{c}(sp));

            % get decay time
            fallSample = dFTrace(spikeLocs{c}(sp):sampleEnd);

            % rise/decay time threshold
            %             [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinder(riseSample, fallSample, exStruct.rate);
            [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinderTriangle(riseSample, fallSample, exStruct.rate);


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

exStruct.cells.zScore = zScore;
exStruct.cells.zScoreThresholded = zScore > zscoreLim;

%% small transient correction

whitenedDF = whiten(exStruct.cells.dF);
riseTimeStructSmall = cell(size(riseTimeStruct));
decayTimeStructSmall = cell(size(riseTimeStruct));
riseTimeIndxStructSmall = cell(size(riseTimeStruct));
decayTimeIndxStructSmall =cell(size(riseTimeStruct));
spikeLocsSmallDF = cell(size(riseTimeStruct));
spikeAmpSmallDF = cell(size(riseTimeStruct));
spikeWidthsSmall = cell(size(riseTimeStruct));

spikeAmpWhitened =[];
for c= 1:exStruct.cellCount
    if exStruct.cells.zScoreThresholded(c) == 1
        dFTraceWhitened = whitenedDF(c,:);
        dFTraceWhitened = movmad(dFTraceWhitened,20);

        dFTrace = exStruct.cells.dF(c,:);

        %% 'spike' detection
        % min peak height
        dF_SD = std(dFTraceWhitened);
        dFMean = mean(dFTraceWhitened);
        minProm = 2 * dF_SD;

        % spike detection on whitened data to get the sample to port over
        % to the raw DF data
        [~,spikeLocsWhitened{c},~] = findpeaks(dFTraceWhitened,"MinPeakHeight",0.4, "MinPeakProminence",minProm);


        % get tau rise time and decay for small spikes in the orignal trace
        if ~isempty(spikeLocsWhitened{c})
            % for each spike
            for sp = 1:length(spikeLocsWhitened{c})

                sampleSt = spikeLocsWhitened{c}(sp)- (exStruct.rate * 2); % take 2 sec before peak onset

                % corection if peak is less than a second before start
                if sampleSt <= 0
                    sampleSt = 1;
                end

                sampleEnd = spikeLocsWhitened{c}(sp)+ (exStruct.rate * 2); % take 2 sec after peak onset

                % corection if sample end is greater than dFtrace end
                if sampleEnd > length(dFTrace)
                    sampleEnd = length(dFTrace);
                end

                % get the peak in the small spike extract
                try
                    dFextract = dFTrace(sampleSt:sampleEnd);
                    [spikeAmpSmallDF{c}(sp), spikeLocsSmallDF{c}(sp), spikeWidthsSmall{c}(sp) ] = findpeaks(dFextract,"MinPeakHeight",0.07,"NPeaks",1);
                    spikeLocsSmallDF{c}(sp) = spikeLocsSmallDF{c}(sp)+sampleSt-1;
                catch
                    spikeAmpSmallDF{c}(sp) = NaN;
                    spikeLocsSmallDF{c}(sp) = NaN;
                    spikeWidthsSmall{c}(sp) = NaN;
                    riseTimeStructSmall{c}(sp) = NaN;
                    decayTimeStructSmall{c}(sp) = NaN;
                    riseTimeIndxStructSmall{c}(sp) = NaN;
                    decayTimeIndxStructSmall{c}(sp) = NaN;
                    break
                end


                % redo the rise and fall times based on 'new' peak
                sampleStrRealigned = (spikeLocsSmallDF{c}(sp))- (exStruct.rate * 1); % take 1 sec before peak onset

                % corection if peak is less than a second before start
                if sampleStrRealigned <= 0
                    sampleStrRealigned = 1;
                end

                sampleEnd = spikeLocsSmallDF{c}(sp)+ (exStruct.rate * 4); % take 4 sec after peak onset

                % corection if sample end is greater than dFtrace end
                if sampleEnd > length(dFTrace)
                    sampleEnd = length(dFTrace);
                end

                % get rise time sample
                riseSample = dFTrace(sampleStrRealigned:spikeLocsSmallDF{c}(sp));

                % get decay time sample
                fallSample = dFTrace(spikeLocsSmallDF{c}(sp):sampleEnd);


                % rise/decay time threshold
                [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinderTriangle(riseSample, fallSample, exStruct.rate);

                % add into structure
                riseTimeStructSmall{c}(sp) = riseTime;
                decayTimeStructSmall{c}(sp) = decayTime;
                riseTimeIndxStructSmall{c}(sp) = spikeLocsSmallDF{c}(sp)-riseTimeIndx;
                decayTimeIndxStructSmall{c}(sp) = spikeLocsSmallDF{c}(sp)+crossIndxFall;
            end
        else
            riseTimeStructSmall{c}(sp) = NaN;
            decayTimeStructSmall{c}(sp) = NaN;
            riseTimeIndxStructSmall{c}(sp) = NaN;
            decayTimeIndxStructSmall{c}(sp) = NaN;
        end
    end

    % remove zeros
    riseTimeStructSmall{c} =  riseTimeStructSmall{c}(riseTimeStructSmall{c}>0);
    decayTimeStructSmall{c} = decayTimeStructSmall{c}(decayTimeStructSmall{c}>0);
    riseTimeIndxStructSmall{c} = riseTimeIndxStructSmall{c}(riseTimeIndxStructSmall{c}>0);
    decayTimeIndxStructSmall{c} = decayTimeIndxStructSmall{c}(decayTimeIndxStructSmall{c}>0);

    % remove nans
    spikeAmpSmallDF{c} =  spikeAmpSmallDF{c}(~isnan(spikeAmpSmallDF{c}));
    spikeLocsSmallDF{c} = spikeLocsSmallDF{c}(~isnan(spikeLocsSmallDF{c}));
    spikeWidthsSmall{c} = spikeWidthsSmall{c}(~isnan(spikeLocsSmallDF{c}));
    riseTimeStructSmall{c} = riseTimeStructSmall{c}(~isnan(riseTimeStructSmall{c}));
    decayTimeStructSmall{c} = decayTimeStructSmall{c}(~isnan(decayTimeStructSmall{c}));
    riseTimeIndxStructSmall{c} = riseTimeIndxStructSmall{c}(~isnan(riseTimeIndxStructSmall{c}));
    decayTimeIndxStructSmall{c} = decayTimeIndxStructSmall{c}(~isnan(decayTimeIndxStructSmall{c}));
end

%% bug fix to remove duplicate spikes, can occur due to whittening and MAD filtering
for vv = 1: length(spikeLocsSmallDF)

    if length(spikeLocsSmallDF{vv}) ~= length(unique(spikeLocsSmallDF{vv}))
        [~, ind] = unique(spikeLocsSmallDF{vv});
        indexToDupes = find(not(ismember(1:numel(spikeLocsSmallDF{vv}),ind)));

        spikeAmpSmallDF{vv}(indexToDupes) = [];
        spikeLocsSmallDF{vv}(indexToDupes) = [];
        spikeWidthsSmall{vv}(indexToDupes) = [];
        riseTimeStructSmall{vv}(indexToDupes) = [];
        decayTimeStructSmall{vv}(indexToDupes) = [];
        riseTimeIndxStructSmall{vv}(indexToDupes) = [];
        decayTimeIndxStructSmall{vv}(indexToDupes) = [];
    end
end
%% remove small spikes too close to large spikes

% for each cell
for c= 1:exStruct.cellCount
    % check if both small and large spikes are not empty
    if ~isempty(spikeLocsSmallDF{c}) && ~isempty(spikeLocs{c})

        % set up temporary data structs to modify
        tempSpikeLocs = spikeLocsSmallDF{c};

        deleteFlag = zeros(1,length(tempSpikeLocs));
        % for each small spike
        for i = 1:length(spikeLocsSmallDF{c})
            % get the frame differences
            frameDiff=  abs(spikeLocs{c}-spikeLocsSmallDF{c}(i));

            % if small spike is close to any large spike then set delete
            % flag

            if any(frameDiff<=frameDiffLimit)
                deleteFlag(i) = 1;
            end
        end

        ind2Del = find(deleteFlag);

        spikeLocsSmallDF{c}(ind2Del) = [];
        riseTimeStructSmall{c}(ind2Del) = [];
        decayTimeStructSmall{c}(ind2Del) = [];
        riseTimeIndxStructSmall{c}(ind2Del) = [];
        decayTimeIndxStructSmall{c}(ind2Del) = [];
        spikeAmpSmallDF{c}(ind2Del) = [];
        spikeWidthsSmall{c}(ind2Del) = [];
    end

    % add firing rate
    if any(riseTimeStructSmall{c})
        firingRateSmall(c) = length(riseTimeStructSmall{c})/ (length(exStruct.xyShifts)/exStruct.rate);
        meanAmpSmall(c) = sum(spikeAmpSmallDF{c})/length(spikeAmpSmallDF{c});
    else
        firingRateSmall(c) = 0;
        meanAmpSmall(c) = 0;
    end

end


%% add everything into exStruct

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

% small spikes
exStruct.spikesSmall.spikeAmp = spikeAmpSmallDF;
exStruct.spikesSmall.spikeLocs = spikeLocsSmallDF;
exStruct.spikesSmall.spikeWidths = spikeWidthsSmall;
exStruct.spikesSmall.riseTime = riseTimeStructSmall;
exStruct.spikesSmall.decayTime = decayTimeStructSmall;
exStruct.spikesSmall.riseTimeIndx = riseTimeIndxStructSmall;
exStruct.spikesSmall.decayTimeIndx = decayTimeIndxStructSmall;
exStruct.spikesSmall.firingRate = firingRateSmall; % firing rate in spk/s
exStruct.spikesSmall.meanSpikeAmp = meanAmpSmall; % mean spike dF amplitude

%% Do the coherence metrics

if sum(exStruct.cells.zScoreThresholded) > 0
    exStruct = waveCellCoherence(exStruct);
end

%% plot everything
zscoreFig = figure('units','normalized','outerposition',[0 0 1 1]);
zscoreAx = gca;
hold on
axis tight
%% plot based on zcore
for c= 1:exStruct.cellCount

    dFTracePlot = exStruct.cells.dF(c,:) + (c-1)*0.4;
    offsetVals(c) = (c-1)*0.4;

    if zScore(c) > zscoreLim
        plot(zscoreAx,dFTracePlot)

        text(zscoreAx,50,dFTracePlot(50),0,num2str(zScore(c)));

        try
            % large spikes
            plot(zscoreAx,spikeLocs{c}(:), dFTracePlot(spikeLocs{c}(:)), 'r*');

            plot(zscoreAx,riseTimeIndxStruct{c} , dFTracePlot(riseTimeIndxStruct{c}), 'b^');
            plot(zscoreAx,decayTimeIndxStruct{c} , dFTracePlot(decayTimeIndxStruct{c}), 'gv');


            % small spikes
            plot(zscoreAx,spikeLocsSmallDF{c}(:), dFTracePlot(spikeLocsSmallDF{c}(:)), '*', 'Color' ,[1 0 1]);

            plot(zscoreAx,riseTimeIndxStructSmall{c} , dFTracePlot(riseTimeIndxStructSmall{c}), 'diamond', 'Color', [0 0 0.5]);
            plot(zscoreAx,decayTimeIndxStructSmall{c} , dFTracePlot(decayTimeIndxStructSmall{c}), 'diamond', 'Color', [0 0.5 0]);


        catch
        end
    else
        plot(zscoreAx,dFTracePlot, '--','Color', [0.5 0.5 0.5]);
        text(zscoreAx,50,dFTracePlot(50),0,num2str(zScore(c)));
    end
end

% plot waves

if sum(exStruct.cells.zScoreThresholded) > 0

    spikesSorted = exStruct.wavesMetrics.spikesSorted;

    for d = 1:max(spikesSorted(:,3))
        wavePoints = spikesSorted(spikesSorted(:,3)==d,:);

        wavePoints(:,2) =  ((wavePoints(:,2))*0.4)-0.4;
        wavePoints(:,3) = [];
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