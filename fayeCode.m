function fayeCode(exStructPath)
% This function completes 'spike' detection on calcium activity of
% idenitified cells. Calculates spike amplitudes, full width at half max,
% tau rise and decays
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs-  exStructPath: filepath for exStruct.mat to process


%% load exStruct
exStruct = load(exStructPath);
exStruct = exStruct.exStruct;

%% analyse each cell

for c= 1:exStruct.cellCount
    dFTrace = exStruct.cells.dF(c,:);

    % min peak height
    dF_SD = std(dFTrace);
    dFMean = mean(dFTrace);
    minHeight = dFMean + (dF_SD *2);
    minProm = 2.5 * dF_SD;

    % spike detection
    [spikeAmp{c},spikeLocs{c},spikeWidths{c}] = findpeaks(dFTrace, "MinPeakProminence",minProm);
    %     [spikeAmp{c},spikeLocs{c},spikeWidths{c}] = findpeaks(dFTrace,"MinPeakHeight",minHeight, "MinPeakProminence",minProm);

    % get tau rise time and decay
    if ~isempty(spikeLocs{c})
        % for eaach spike
        for sp = 1:length(spikeLocs{c})

            sampleSt = spikeLocs{c}(sp)- (exStruct.rate * 4); % take 4 sec before peak onset

            % corection if peak is less than a second before start
            if sampleSt < 0
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
            [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinder(riseSample, fallSample, exStruct.rate);

            %             plot(dFTrace)
            %             hold on
            %             plot(spikeLocs{c}(sp), dFTrace(spikeLocs{c}(sp)), '*');
            %
            %             plot(spikeLocs{c}(sp)-riseTimeIndx , dFTrace(spikeLocs{c}(sp) -riseTimeIndx), 'b*');
            %             plot(spikeLocs{c}(sp)+crossIndxFall , dFTrace(spikeLocs{c}(sp) +crossIndxFall), 'g*');
            %
            %             hold off

            % add into structure
            riseTimeStruct{c}(sp) = riseTime;
            decayTimeStruct{c}(sp) = decayTime;
            riseTimeIndxStruct{c}(sp) = spikeLocs{c}(sp)-riseTimeIndx;
            decayTimeIndxStruct{c}(sp) = spikeLocs{c}(sp)+crossIndxFall;
        end
    else
    end
end

%% add everything into exStruct

exStruct.spikes.spikeAmp = spikeAmp;
exStruct.spikes.spikeLocs = spikeLocs;
exStruct.spikes.spikeWidths = spikeWidths;
exStruct.spikes.riseTime = riseTimeStruct;
exStruct.spikes.decayTime = decayTimeStruct;
exStruct.spikes.riseTimeIndx = riseTimeIndxStruct;
exStruct.spikes.decayTimeIndx = decayTimeIndxStruct;

% save data
save(exStructPath, "exStruct", '-v7.3');
end