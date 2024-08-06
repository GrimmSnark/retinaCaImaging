function concatSpikeDataFaye()
% pulls all the spike data into a single excel sheet

%  [path] = uigetdir('', 'Pick a Directory');
[path] = '\\campus\rdw\ion10\10\retina\data\Savage\calcium\Faye\July 2024 files';


exStructPaths = dir([path '\**\*exStruct.mat']);

disp(['Found ' num2str(length(exStructPaths)) ' experiment folders, moving on....'])

firingRates = [];
meanSpikeAmp = [];
meanExStructPath = [];
cellNoMean = [];
exStructTablePath = [];

firingRatesSmallSpikes = [];
meanSpikeAmpSmallSpikes = [];
meanExStructPathSmallSpikes = [];
cellNoMeanSmallSpikes = [];
exStructTablePathSmallSpikes = [];

cellCoherencePaths =[];
cellCoherenceCellNo =[];
cellCoherence =[];

waveCoherencePaths = [];
waveCoherenceWaveNo = [];
waveCoherence = [];

%% spike metrics
spikeAmp = [];
spikeWidth = [];
riseTime = [];
decayTime = [];
cellNo = [];
spikeNo = [];

spikeSmallAmp = [];
spikeSmallWidth = [];
riseTimeSmall = [];
decayTimeSmall = [];
cellNoSmall = [];
spikeSmallNo = [];

% for all exStructs
for i = 1:length(exStructPaths)
    tempTable = table;
    currentExPath = fullfile(exStructPaths(i).folder, exStructPaths(i).name);

    load(currentExPath);

    if isfield(exStruct, 'spikes')
        for c = 1:exStruct.cellCount

            % large spike stuff
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


            % small spike stuff
            currAmpsSmall = exStruct.spikesSmall.spikeAmp{c}';
            currWidthsSmall = exStruct.spikesSmall.spikeWidths{c}';
            currRiseSmall = exStruct.spikesSmall.riseTime{c}';
            currDecaySmall = exStruct.spikesSmall.decayTime{c}';
            currCellSmall = repmat(c,1,length(currDecaySmall))';

            spikeSmallAmp = [spikeSmallAmp ;currAmpsSmall];
            spikeSmallWidth = [ spikeSmallWidth; currWidthsSmall];
            riseTimeSmall = [riseTimeSmall; currRiseSmall];
            decayTimeSmall = [decayTimeSmall ; currDecaySmall];
            cellNoSmall = [cellNoSmall ;currCellSmall];
            spikeSmallNo = [spikeSmallNo; (1:length(currCellSmall))'];
        end

        exStructTablePath = [exStructTablePath ;repmat({currentExPath}, length([exStruct.spikes.spikeAmp{:}]),1)];
        exStructTablePathSmallSpikes = [exStructTablePathSmallSpikes ;repmat({currentExPath}, length([exStruct.spikesSmall.spikeAmp{:}]),1)];

        %% mean metrics
        meanExStructPath = [meanExStructPath ;repmat({currentExPath}, exStruct.cellCount,1)];
        cellNoMean = [cellNoMean; (1: exStruct.cellCount)'];
        firingRates = [firingRates; exStruct.spikes.firingRate' ];
        meanSpikeAmp = [meanSpikeAmp; exStruct.spikes.meanSpikeAmp'];

        meanExStructPathSmallSpikes = [meanExStructPathSmallSpikes ;repmat({currentExPath}, exStruct.cellCount,1)];
        cellNoMeanSmallSpikes = [cellNoMeanSmallSpikes; (1: exStruct.cellCount)'];
        firingRatesSmallSpikes = [firingRatesSmallSpikes; exStruct.spikesSmall.firingRate' ];
        meanSpikeAmpSmallSpikes = [meanSpikeAmpSmallSpikes; exStruct.spikesSmall.meanSpikeAmp'];

        %% cell coherence metrics
        cellCoherencePaths = [cellCoherencePaths ;repmat({currentExPath}, length(exStruct.wavesMetrics.cellCoherence) ,1)];
        cellCoherenceCellNo = [cellCoherenceCellNo; exStruct.wavesMetrics.cellCoherence(:,1)];
        cellCoherence = [cellCoherence; exStruct.wavesMetrics.cellCoherence(:,2)];

        %% wave coherence metrics
        waveCoherencePaths = [waveCoherencePaths ; repmat({currentExPath}, length(exStruct.wavesMetrics.waveCoherence) ,1)];
        waveCoherenceWaveNo = [ waveCoherenceWaveNo ; (1:length(exStruct.wavesMetrics.waveCoherence))'];
        waveCoherence = [waveCoherence; exStruct.wavesMetrics.waveCoherence'];

    else

        disp([currentExPath ' does not contain spikes field....moving on']);
    end

end

grandTable = table(exStructTablePath, cellNo,  spikeNo, spikeAmp, spikeWidth, riseTime, decayTime);
grandTableSmall = table(exStructTablePathSmallSpikes, cellNoSmall,  spikeSmallNo, spikeSmallAmp, spikeSmallWidth, riseTimeSmall, decayTimeSmall);

meanTable = table(meanExStructPath, cellNoMean, firingRates, meanSpikeAmp);
meanTableSmall = table(meanExStructPathSmallSpikes, cellNoMeanSmallSpikes, firingRatesSmallSpikes, meanSpikeAmpSmallSpikes);

cellCoherenceTable = table(cellCoherencePaths, cellCoherenceCellNo, cellCoherence);
waveCoherenceTable = table(waveCoherencePaths, waveCoherenceWaveNo, waveCoherence);

% clean tables
meanTable(meanTable.meanSpikeAmp == 0,:)=[];
meanTableSmall(meanTableSmall.meanSpikeAmpSmallSpikes == 0,:)=[];


% saving

writetable(grandTable, fullfile(path, '\summaryExcel.xlsx'),'Sheet','Big Spike Grand Table');
writetable(meanTable, fullfile(path, '\summaryExcel.xlsx'),'Sheet','Big Spike Mean Table');
writetable(grandTableSmall, fullfile(path, '\summaryExcel.xlsx'),'Sheet','Small Spike Grand Table');
writetable(meanTableSmall, fullfile(path, '\summaryExcel.xlsx'),'Sheet','Small Spike Mean Table');
writetable(cellCoherenceTable, fullfile(path, '\summaryExcel.xlsx'),'Sheet','Cell Coherence Table');
writetable(waveCoherenceTable, fullfile(path, '\summaryExcel.xlsx'),'Sheet','Wave Coherence Table');

end