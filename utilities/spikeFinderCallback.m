function app = spikeFinderCallback(app)

%% dialog set up
prompt = {'Rise Time Limit (sec):','Decay Time Limit (sec):'};
dlgtitle = 'Input rise and decay time limits';
dims = [1 35];
definput = {'8','8'};
dlgResponse = inputdlg(prompt,dlgtitle,dims,definput);


%% do the recalculation
c = app.cellNum;
dFTrace = app.exStruct.cells.dF(c,:);

% min peak height
dF_SD = std(dFTrace);
dFMean = mean(dFTrace);
minHeight = dFMean + (dF_SD *2.5);
minProm = 2.5 * dF_SD;

% spike detection
[spikeAmp,spikeLocs,spikeWidths] = findpeaks(dFTrace, "MinPeakProminence",minProm);

% get tau rise time and decay
if ~isempty(spikeLocs)
    % for eaach spike
    for sp = 1:length(spikeLocs)

        sampleSt = spikeLocs(sp)- (app.exStruct.rate * str2double(dlgResponse{1})); 

        % corection if peak is less than a second before start
        if sampleSt <= 0
            sampleSt = 1;
        end

        sampleEnd = spikeLocs(sp)+ (app.exStruct.rate * str2double(dlgResponse{2}));

        % corection if sample end is greater than dFtrace end
        if sampleEnd > length(dFTrace)
            sampleEnd = length(dFTrace);
        end

        % get rise time
        riseSample = dFTrace(sampleSt:spikeLocs(sp));

        % get decay time
        fallSample = dFTrace(spikeLocs(sp):sampleEnd);

        % rise/decay time threshold
        [riseTime, decayTime, riseTimeIndx, crossIndxFall] = riseTimeDecayFinderTriangle(riseSample, fallSample, app.exStruct.rate);


        % add into structure
        riseTimeStruct(sp) = riseTime;
        decayTimeStruct(sp) = decayTime;
        riseTimeIndxStruct(sp) = spikeLocs(sp)-riseTimeIndx;
        decayTimeIndxStruct(sp) = spikeLocs(sp)+crossIndxFall;
    end
else
end

% add firing rate
app.exStruct.spikes.firingRate(c) = length(riseTimeStruct)/ (length(app.exStruct.xyShifts)/app.exStruct.rate);
app.exStruct.spikes.meanSpikeAmp.meanAmp(c) = sum(spikeAmp)/length(spikeAmp);

%% overwrite app structure
app.exStruct.spikes.spikeAmp{c} = spikeAmp;
app.exStruct.spikes.spikeLocs{c} = spikeLocs;
app.exStruct.spikes.spikeWidths{c} = spikeWidths;
app.exStruct.spikes.riseTime{c} = riseTimeStruct;
app.exStruct.spikes.decayTime{c} = decayTimeStruct;
app.exStruct.spikes.riseTimeIndx{c} = riseTimeIndxStruct;
app.exStruct.spikes.decayTimeIndx{c} = decayTimeIndxStruct;

exStruct = app.exStruct;

% save data
exStructPath = [app.exStruct.filePath(1:end-5) '.nd2'];
save(exStructPath, "exStruct", '-v7.3');


end