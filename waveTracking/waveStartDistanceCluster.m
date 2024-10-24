function exStruct = waveStartDistanceCluster(exStruct)

clusterInfo = exStruct.clusterInfo;
clusterVerts = clusterInfo.clusterPolys;
waveMetrics = exStruct.waves.relDistanceTab;

clusterMap = zeros(size(exStruct.dF,[1 2]));
% build cluster image

for i = 1:length(clusterVerts)
    currentCluster = clusterVerts{i};

    % image at full resolution
    tempIm = poly2mask(currentCluster(:,1), currentCluster(:,2), exStruct.image.pixelNum, exStruct.image.pixelNum);

    % downsample to dF resolution
    downSampleFactor = size(exStruct.dF,1)  /  exStruct.image.pixelNum;

    tempImResize = imresize(tempIm,downSampleFactor, "nearest");
    clusterMap = clusterMap + tempImResize;
end

[clusterPixelsX, clusterPixelsY] = find(clusterMap);
% run through waves and get the distance to the closest cluster pixel

for w = 1:length(exStruct.waves.centerPerFrame)
    [closestClusterPerWaveStart(w), minIn(w)]= min(pdist2(exStruct.waves.centerPerFrame{w}(1,:), [clusterPixelsX, clusterPixelsY],"euclidean"));

end

exStruct.waves.relDistanceTab.closestClusterPixelWaveStart = closestClusterPerWaveStart';
end