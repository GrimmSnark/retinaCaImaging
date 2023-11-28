function concatSpikeDataFaye()
% pulls all the spike data into a single excel sheet

%  [path] = uigetdir('', 'Pick a Directory');
[path] = 'I:\Faye';


exStructPaths = dir([path '\**\*exStruct.mat']);

disp(['Found ' num2str(length(exStructPaths)) ' experiment folders, moving on....'])

firingRates = [];
meanSpikeAmp = [];
meanExStructPath = [];
cellNoMean = [];
exStructTablePath = [];

%% spike metrics
spikeAmp = [];
spikeWidth = [];
riseTime = [];
decayTime = [];
cellNo = [];
spikeNo = [];

% for all exStructs
for i = 1:length(exStructPaths)
    tempTable = table;
    currentExPath = fullfile(exStructPaths(i).folder, exStructPaths(i).name);

    load(currentExPath);

    if isfield(exStruct, 'spikes')
        for c = 1:exStruct.cellCount
            currAmps = exStruct.spikes.spikeAmp{c}';
            currWidths = exStruct.spikes.spikeWidths{c}';
            currRise = exStruct.spikes.riseTime{c}';
            currDecay = exStruct.spikes.decayTime{c}';
            currCell = repmat(c,1,length(currDecay))';

            spikeAmp = [spikeAmp ;currAmps];
            spikeWidth = [ spikeWidth; currWidths];
            riseTime = [riseTime; currRise];
            decayTime = [decayTime ; currDecay];
            cellNo = [cellNo ;currCell];
            spikeNo = [spikeNo; (1:length(currCell))'];

        end

        exStructTablePath = [exStructTablePath ;repmat({currentExPath}, length([exStruct.spikes.spikeAmp{:}]),1)];

        %% mean metrics
        meanExStructPath = [meanExStructPath ;repmat({currentExPath}, exStruct.cellCount,1)];
        cellNoMean = [cellNoMean; (1: exStruct.cellCount)'];
        firingRates = [firingRates; exStruct.spikes.firingRate' ];
        meanSpikeAmp = [meanSpikeAmp; exStruct.spikes.meanSpikeAmp'];

    else

        disp([currentExPath 'does not contain spikes field....moving on']);
    end

end

grandTable = table(exStructTablePath, cellNo,  spikeNo, spikeAmp, spikeWidth, riseTime, decayTime);
meanTable = table(meanExStructPath, cellNoMean, firingRates, meanSpikeAmp);

% saving

writetable(grandTable, fullfile(path, '\summaryExcel.xlsx'),'Sheet',1);
writetable(meanTable, fullfile(path, '\summaryExcel.xlsx'),'Sheet',2);

end