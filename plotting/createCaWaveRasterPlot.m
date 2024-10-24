function createCaWaveRasterPlot(exStruct)

scaleFactor = 1/4;
%% downsample
% dwnSampleDF = imresize(exStruct.dF, 0.125);
% dwnSampleBoundary = imresize(exStruct.waves.retinaBoundMask, 0.125);
dwnSampleDF = imresize(exStruct.dF, scaleFactor);
dwnSampleBoundary = imresize(exStruct.waves.retinaBoundMask, scaleFactor);

linDF = [];
linDF = gpuArray(reshape(dwnSampleDF, [], size(dwnSampleDF,3)));
linBoundary = gpuArray(reshape(dwnSampleBoundary,[],1));
%% run through pixels for spike times

rasterIm= [];
parfor i = 1: length(linDF)
    intensityThres = mean(linDF(i,:)) + (std(linDF(i,:)) *3);

    rasterIm(i,:) = linDF(i,:)>intensityThres;

end

%% align each spike to the waveID

waveTab = exStruct.waves.waveTable(:,[4 5 8]);

% go through frames
for fr = 1:size(rasterIm,2)
    pixelSpksPerFrame = find(rasterIm(:,fr) == 1); % find index of 'spikes'
    tempIm = zeros(size(dwnSampleDF,[1 2]));
    tempIm(pixelSpksPerFrame) = 1;
    % get the available 

    % for each pixel
    for pix = 1:size()

    end
end

rasterIm = rasterIm(logical(linBoundary),:);
rasterIm = logical(rasterIm);
plotCaRaster(rasterIm);
plotSpikeRaster(rasterIm, 'PlotType',  'vertline');
% imshow(rasterIm);
end