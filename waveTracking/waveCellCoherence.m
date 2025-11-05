function exStruct = waveCellCoherence(exStruct)

%% defaults

waveFrameLim = 5;
spikeNo4Waves = 2;

%%
% unwrap all spikes and match up cell IDs
spikes = exStruct.spikes;

spikeUnwrapped = [];
for i =1:length(spikes.spikeLocs)
    if exStruct.cells.zScoreThresholded(i) == 1
        tempSpikes = spikes.spikeLocs{i}';
        tempSpikes(:,2) = ones(length(tempSpikes),1)*i;
        spikeUnwrapped = [spikeUnwrapped;tempSpikes];
    end
end

% sort spikes by frames
spikesSorted = sortrows(spikeUnwrapped,1);
spikesSortedDiff = spikesSorted(2:end,1)-spikesSorted(1:end-1,1);

% find wave breaks
waveEnds = find(spikesSortedDiff>waveFrameLim);
waveEnds(end+1,1) = length(spikesSorted);
waveStarts = [1 ; waveEnds+1];
waveStarts = waveStarts(1:end-1);

% align waves to spikes
waveIDs = [];
for i = 1:length(waveStarts)
    waveIDs(waveStarts(i):waveEnds(i)) = i;
end

% clean waves to remove 'single spikes'
for x= 1:max(waveIDs)
    if sum(waveIDs == x) <= spikeNo4Waves
        waveIDs(waveIDs == x) = NaN;
    end
end

uniqueNums = rmmissing(unique(waveIDs));
waveIDs = changem(waveIDs,1:length(uniqueNums),uniqueNums);
waveIDs(isnan(waveIDs))= 0;

spikesSorted(:,3) = waveIDs;

%% calculate wave "coherence"

if sum(~isnan(waveIDs)) > 0
    for i = 1:max(waveIDs)
        waveCoh(i) = sum(spikesSorted(:,3)==i)/length(unique(spikesSorted(:,2)));
    end
else
    waveCoh = [];
end

%% cell "coherence"

count = 1;
cellNums = unique(spikesSorted(:,2));
% for all cells which spike
for i = 1:length(cellNums)
    tempCellNum = cellNums(i);
    tempCellSpikes = [];

    if exStruct.cells.zScoreThresholded(tempCellNum) == 1
        tempCellSpikes = spikesSorted(:,2)==tempCellNum;
        tempSpikeList = spikesSorted(tempCellSpikes,:);
        cellCoherence(count,1)= tempCellNum;
        cellCoherence(count,2)= 1 - (sum(~(tempSpikeList(:,3)))/length(tempSpikeList));
        count = count+1;
    end
end

%% put everything in structure
exStruct.wavesMetrics.spikesSorted = spikesSorted;
exStruct.wavesMetrics.waveCoherence = waveCoh;
exStruct.wavesMetrics.cellCoherence = cellCoherence;


end