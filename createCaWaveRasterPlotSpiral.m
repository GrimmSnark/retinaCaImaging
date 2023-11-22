function createCaWaveRasterPlotSpiral(exStructPath, saveType, downsampleFactor)

if nargin < 1 || isempty(exStructPath)
    [file, path] = uigetfile({'*.mat'},...
        'Experiment Struct Selector');

    exStructPath = fullfile(path,file);
end

if nargin <2 || isempty(saveType)
    saveType = '.tif';
end

if nargin < 3 || isempty(downsampleFactor)
    downsampleFactor = 1;
end

%% load in exStruct
exStruct = load(exStructPath);
exStruct = exStruct.exStruct;
disp('Loaded in exStruct.mat')

%% load in waveCol movie

[saveFolder, name] = fileparts(exStructPath);
waveColPath = fullfile(saveFolder, [name(1:end-9) '_waveCol.tif']);
waveIm = read_Tiffs(waveColPath,[],1);

waveIm = imresize(waveIm,downsampleFactor,"nearest");

%% contruct filtering cmap
waveCmap(1,:) = [ 0 0 0];
waveCmap(2,:) = [1 1 1];
waveCols = distinguishable_colors(length(unique(exStruct.waves.waveTable.waveNumber)),{'w','k'});

waveCmap = [waveCmap ; waveCols];

%% convert rgb to index image
indIm = zeros(size(waveIm,1), size(waveIm,2), size(waveIm,4));

for i = 1:size(waveIm,4)
    indIm(:,:,i) = rgb2ind(waveIm(:,:,:,i),waveCmap);
end

indIm(indIm<2) = 0;
indIm(indIm>=2) = indIm(indIm>=2)-1;

indImReshape = reshape(indIm, [], size(indIm,3));

%% get retinal center location and spiral indexing

retinaBoundMask = imresize(exStruct.waves.retinaBoundMask,downsampleFactor,"nearest");

centerCoor = regionprops(retinaBoundMask, 'Centroid');
centerCoor = round(centerCoor.Centroid);
spiralIndexing = spiral_generic(size(waveIm, [1 2]), centerCoor);
spiralIndexLinear = reshape(spiralIndexing, [],1);

indImReshape = indImReshape(spiralIndexLinear,:);


%% restrict to pixels inside boundary
dwnSampleBoundary = imresize(exStruct.waves.retinaBoundMask, downsampleFactor, "nearest");
linBoundary = reshape(dwnSampleBoundary,[],1);

% reorder based on spiral
linBoundary = linBoundary(spiralIndexLinear);

indImReshape = indImReshape(logical(linBoundary),:);

%% plot pixel x timeseries image
% waveMapDisplay = [1 1 1 ; waveCols];
% 
% figH = figure('units','normalized','outerposition',[0 0 1 1]);
% imagesc(indImReshape);
% colormap(waveMapDisplay);
% xlim([0 size(indImReshape,2)]);
% ylim([0 size(indImReshape,1)]);
% 
% xlabel('Time in seconds');
% ylabel(['Pixel ROI no. (Resolution ' num2str(round(exStruct.downsampledRes)) 'um x ' num2str(round(exStruct.downsampledRes))  'um)']);
% 
% 
% title(name,Interpreter="none");
% tightfig;

figH = figure('units','normalized','outerposition',[0 0 1 1]);
hold on

for i = 1:max(exStruct.waves.waveTable.waveNumber)
    spikes = indImReshape ==i;
 [xPoints, yPoints] = plotCaRaster(spikes, waveCols(i,:));
end

xlim([0 size(indImReshape,2)]);
ylim([0 size(indImReshape,1)]);

xlabel('Time in seconds');
ylabel(['Pixel ROI no. (Resolution ' num2str(round(exStruct.downsampledRes)/downsampleFactor) 'um x ' num2str(round(exStruct.downsampledRes)/downsampleFactor)  'um)']);


title(name,Interpreter="none");
tightfig;

%% save

if ~strcmp(saveType,'.eps')
    saveas(figH, fullfile(saveFolder, [name(1:end-9) '_waveRasterSpiral' saveType]));
else
    exportgraphics(gca, fullfile(saveFolder, [name(1:end-9) '_waveRasterSpiral' saveType]), "ContentType","vector");
end

close;
end