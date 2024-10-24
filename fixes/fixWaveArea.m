function exStruct = fixWaveArea(exStruct)


tempArea = [];
tempAreaMicron =[];
for w = 1:max(exStruct.waves.waveTable.waveNumber)

    % get all objects for a single waves
    subTabIndx = exStruct.waves.waveTable.waveNumber == w;
    subTable = exStruct.waves.waveTable(subTabIndx,:);

    % get all the pixels involved in wave
    wavePixels = unique(cat(1,subTable.PixelIdxList{:}));
    [wavePixX, wavePixY] = ind2sub([512 512], wavePixels);

    % get wave extent in pixels
    tempArea(w) = length(wavePixels);
    tempAreaMicron(w) = length(wavePixels) *  exStruct.downsampledRes;
end
exStruct.waves.waveArea = tempArea;
exStruct.waves.waveAreaMicron = tempAreaMicron;

end