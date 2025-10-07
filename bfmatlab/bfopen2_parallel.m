function [result] = bfopen2_parallel(id, varargin)
% bfopen2_parallel - Parallel version of bfopen2 using Bio-Formats
% 
% This version parallelizes the loading of image planes using parfor.
% Make sure to initialize MATLAB's parallel pool before running.
% Requires Bio-Formats installed and configured.

autoloadBioFormats = 1;
stitchFiles = 0;

% -- Bio-Formats initialization --
status = bfCheckJavaPath(autoloadBioFormats);
assert(status, 'Bio-Formats library is not on the path.');

if nargin == 0 || exist(id, 'file') == 0
    [file, path] = uigetfile(bfGetFileExtensions, 'Choose a file to open');
    id = [path file];
    if isequal(path, 0) || isequal(file, 0), return; end
end

bfInitLogging();

% Initialize memoized reader
r = bfGetReader();
r = loci.formats.Memoizer(r);
r.setId(id);
r.close();
r.setId(id);

if nargin >=4
    planeSize = javaMethod('getPlaneSize', 'loci.formats.FormatTools', ...
                           r, varargin{3}, varargin{4});
else
    planeSize = javaMethod('getPlaneSize', 'loci.formats.FormatTools', r);
end

if planeSize / (1024)^3 >= 2
    error('Image plane too large. Try reading in tiles.');
end

numSeries = r.getSeriesCount();
result = cell(numSeries, 2);

globalMetadata = r.getGlobalMetadata();

for s = 1:numSeries
    fprintf('Reading series #%d...\n', s);
    r.setSeries(s - 1);
    pixelType = r.getPixelType();
    bpp = javaMethod('getBytesPerPixel', 'loci.formats.FormatTools', pixelType);
    bppMax = power(2, bpp * 8);
    numImages = r.getImageCount();
    sizeT = r.getSizeT();

    % -- Parallel Processing Chunk Setup --
    nWorkers = min(parcluster('local').NumWorkers, numImages);
    chunkSize = ceil(numImages / nWorkers);

    imageChunks = cell(1, nWorkers);
    colorMapChunks = cell(1, nWorkers);

    parfor w = 1:nWorkers
        localReader = bfGetReader();
        localReader = loci.formats.Memoizer(localReader);
        localReader.setId(id);
        localReader.setSeries(s - 1);

        tStart = (w - 1) * chunkSize + 1;
        tEnd   = min(w * chunkSize, numImages);
        localImageList = cell(tEnd - tStart + 1, 2);
        localColorMaps = cell(tEnd - tStart + 1, 1);

        for i = tStart:tEnd
            arr = bfGetPlane(localReader, i, varargin{:});
            
            % Color map
            if bpp == 1
                cmap = localReader.get8BitLookupTable()';
            else
                cmap = localReader.get16BitLookupTable()';
            end
            if ~isempty(cmap)
                cmap = single(cmap);
                cmap(cmap < 0) = cmap(cmap < 0) + bppMax;
                cmap = cmap / (bppMax - 1);
            end

            % Label
            label = id;
            if numImages > 1
                zct = localReader.getZCTCoords(i - 1);
                qt = int2str(zct(3) + 1);
                label = sprintf('%s; T=%s/%d', label, qt, sizeT);
            end

            localImageList{i - tStart + 1, 1} = arr;
            localImageList{i - tStart + 1, 2} = label;
            localColorMaps{i - tStart + 1} = cmap;
        end

        localReader.close();
        imageChunks{w} = localImageList;
        colorMapChunks{w} = localColorMaps;
    end

    % Combine chunks into single output
    imageList = vertcat(imageChunks{:});
    colorMaps = vertcat(colorMapChunks{:});

    % Save results
    result{s, 1} = imageList;

    % Series metadata
    seriesMetadata = r.getSeriesMetadata();
    javaMethod('merge', 'loci.formats.MetadataTools', ...
               globalMetadata, seriesMetadata, 'Global ');
    result{s, 2} = seriesMetadata;
    result{s, 3} = colorMaps;
    result{s, 4} = r.getMetadataStore();

    fprintf('Done reading series #%d\n', s);
end

r.close();

end