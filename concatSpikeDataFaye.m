function concatSpikeDataFaye(exStructFolder)
% pulls all the spike data into a single excel sheet

if nargin < 1 || isempty(exStructFolder)

 [exStructFolder] = uigetdir('', 'Pick a Directory');
% [path] = '\\campus\rdw\ion10\10\retina\data\Savage\calcium\Faye\July 2024 files';
end


exStructPaths = dir([exStructFolder '\**\*exStruct.mat']);

disp(['Found ' num2str(length(exStructPaths)) ' experiment folders, moving on....'])

firingRates = [];
meanSpikeAmp = [];
meanExStructPath = [];
cellNoMean = [];
exStructTablePath = [];

cellCoherencePaths =[];
cellCoherenceCellNo =[];
cellCoherence =[];

waveCoherencePaths = [];
waveCoherenceWaveNo = [];
waveCoherence = [];

zScoresPaths = [];
zScoreCells = [];
zScores = [];

%% spike metrics
spikeAmp = [];
spikeWidth = [];
riseTime = [];
decayTime = [];
cellNo = [];
spikeNo = [];
waveInclusion = [];


% for all exStructs
for i = 1:length(exStructPaths)
    tempTable = table;
    currentExPath = fullfile(exStructPaths(i).folder, exStructPaths(i).name);

    load(currentExPath);
    zInclude = exStruct.cells.zScoreThresholded; % filter by z score 

    if isfield(exStruct, 'spikes')
        for c = 1:exStruct.cellCount

            if exStruct.cells.zScoreThresholded(c) == 1
                % large spike stuff
                currAmps = exStruct.spikes.spikeAmp{c}';
                currWidths = exStruct.spikes.spikeWidths{c}';
                currRise = exStruct.spikes.riseTime{c}';
                currDecay = exStruct.spikes.decayTime{c}';
                currCell = repmat(c,1,length(currDecay))';
                waveInclusionTemp = exStruct.wavesMetrics.spikesSorted(exStruct.wavesMetrics.spikesSorted(:,2)==c,3); % wave inclusion per spike
            
                spikeAmp = [spikeAmp ;currAmps];
                spikeWidth = [ spikeWidth; currWidths];
                riseTime = [riseTime; currRise];
                decayTime = [decayTime ; currDecay];
                cellNo = [cellNo ;currCell];
                spikeNo = [spikeNo; (1:length(currCell))'];
                waveInclusion = [waveInclusion; waveInclusionTemp];
            end
        end

        tempSpikeNum = length([exStruct.spikes.spikeAmp{exStruct.cells.zScoreThresholded}]); % gives number of spikes from zScore thresholded cells
        exStructTablePath = [exStructTablePath ;repmat({currentExPath}, tempSpikeNum,1)];
        %% mean metrics
        meanExStructPath = [meanExStructPath ;repmat({currentExPath}, sum(zInclude),1)];
        cellNoMean = [cellNoMean; find(exStruct.cells.zScoreThresholded)];
        firingRates = [firingRates; exStruct.spikes.firingRate(zInclude)' ];
        meanSpikeAmp = [meanSpikeAmp; exStruct.spikes.meanSpikeAmp(zInclude)'];

        %% cell coherence metrics
        cellCoherencePaths = [cellCoherencePaths ;repmat({currentExPath}, sum(zInclude) ,1)];
        cellCoherenceCellNo = [cellCoherenceCellNo; exStruct.wavesMetrics.cellCoherence(:,1)];
        cellCoherence = [cellCoherence; exStruct.wavesMetrics.cellCoherence(:,2)];

        %% wave coherence metrics
        waveCoherencePaths = [waveCoherencePaths ; repmat({currentExPath}, length(exStruct.wavesMetrics.waveCoherence) ,1)];
        waveCoherenceWaveNo = [ waveCoherenceWaveNo ; (1:length(exStruct.wavesMetrics.waveCoherence))'];
        waveCoherence = [waveCoherence; exStruct.wavesMetrics.waveCoherence'];


        zScoresPaths = [zScoresPaths ; repmat({currentExPath}, length(exStruct.cells.zScore) ,1)];
        zScoreCells = [zScoreCells ; (1:exStruct.cellCount)'];
        zScores = [zScores; exStruct.cells.zScore];
    else

        disp([currentExPath ' does not contain spikes field....moving on']);
    end

end

grandTable = table(exStructTablePath, cellNo,  spikeNo, spikeAmp, spikeWidth, riseTime, decayTime, waveInclusion);
grandTableSyncedSpikes = grandTable(grandTable.waveInclusion >0,:);
grandTableAsynchedSpikes = grandTable(grandTable.waveInclusion ==0,:);

meanTable = table(meanExStructPath, cellNoMean, firingRates, meanSpikeAmp);

cellCoherenceTable = table(cellCoherencePaths, cellCoherenceCellNo, cellCoherence);
waveCoherenceTable = table(waveCoherencePaths, waveCoherenceWaveNo, waveCoherence);


zScoreTable = table(zScoresPaths,zScoreCells,zScores);

% clean tables
meanTable(meanTable.meanSpikeAmp == 0,:)=[];

% saving

writetable(grandTable, fullfile(exStructFolder, '\summaryExcel.xlsx'),'Sheet','Spike Grand Table');
writetable(grandTableSyncedSpikes,fullfile(exStructFolder, '\summaryExcel.xlsx'),'Sheet','Synced Spike Grand Table');
writetable(grandTableAsynchedSpikes, fullfile(exStructFolder, '\summaryExcel.xlsx'),'Sheet','Asynced Spike Grand Table');
writetable(meanTable, fullfile(exStructFolder, '\summaryExcel.xlsx'),'Sheet','Spike Mean Table');
writetable(zScoreTable, fullfile(exStructFolder, '\summaryExcel.xlsx'),'Sheet','Cell ZScores');
writetable(cellCoherenceTable, fullfile(exStructFolder, '\summaryExcel.xlsx'),'Sheet','Cell Coherence Table');
writetable(waveCoherenceTable, fullfile(exStructFolder, '\summaryExcel.xlsx'),'Sheet','Wave Coherence Table');

end