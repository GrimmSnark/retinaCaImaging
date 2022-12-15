function trackCalciumWaves(filepathDF, thresholdVal, waveMinSize)
% Tracks calcium waves and summarises frequency, speed , trajectory, wave
% size etc. Requires dF/F pixelwise movie to have been created from
% prepRetinaCalcium
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Input- filepathDF: filepath for dF/F tif stack
%
%        thresholdVal: threshold value for blob detection (0 - 1 range)
%                      OPTIONAL, DEFAULT = 0.01
%
%        waveMinSize: Calcium wave minimum size in pixel area
%                     OPTIONAL, DEFAULT = 300

%% defaults

if nargin < 2 || isempty(thresholdVal)
   thresholdVal = 0.01; 
end

if nargin < 3 || isempty(waveMinSize)
    waveMinSize = 300;
end

%% load in dF movie
tifStack = read_Tiffs(filepathDF);
[folderPath, name] = fileparts(filepathDF);

%% load in exStruct
exStruct = load(fullfile(folderPath,[name(1:end-5) '_ExStruct.mat']));

try
    exStruct = exStruct.exStruct;
catch
    exStruct = exStruct.metaData;
end

%% get the downsampled res factor
downSampleFactor = exStruct.image.pixelNum/size(tifStack,1);
exStruct.downsampledRes = downSampleFactor*exStruct.image.pixelSize;

      
%% gaus blur
tifGaus = imgaussfilt(tifStack,4);

%% remove background
percent2Subtract = 5;
frameBrightness = squeeze(mean(tifGaus, [1 2]));
highpassFilteredTrace = baselinePercentileFilter(frameBrightness,3,10,percent2Subtract);

% plot(frameBrightness);
% hold on
% plot(highpassFilteredTrace);

for x = 1:size(tifGaus,3)
    tifGauSub(:,:,x) = tifGaus(:,:,x) - uint16(highpassFilteredTrace(x));
end

%% move to FIJI to use their threshold
stackImp = MIJ.createImage('tifs', tifGauSub, 1);

% get frame brightness and the 3/4 brightest frame number
frameBrightness = squeeze(mean(tifGauSub,[1 2]));

[frameVal, sortedIndx ]= sort(frameBrightness, 'ascend');
frame2Set = round(length(frameBrightness) * 0.90);

MIJ.setSlice(sortedIndx(frame2Set));
MIJ.run('Threshold...','setAutoThreshold=Mean');
MIJ.run("Convert to Mask", "method=Default background=Dark black");
MIJ.run("Fill Holes", "stack");

tifTreshFilled = im2uint8(rescale(MIJ.getImage('tifs')));

% Clean up windows
stackImp.changes = false;
stackImp.close;
MIJ.closeAllWindows;

threshBinary = imbinarize(tifTreshFilled);
disp('Thresholding done...');

% %% threshold
% tifThresh = imbinarize(tifGauSub, thresholdVal);
% tifThresh = im2uint8(tifThresh);

%% fill holes to make wave detection easier
% for cc = 1:size(tifThresh,3)
%     tifTreshFilled(:,:,cc) = imfill(tifThresh(:,:,cc));
% end

disp('Holes filled');

% threshBinary = imbinarize(tifTreshFilled);

%%
waveTable = [];
% opticFlow =opticalFlowHS;
for fr = 1:size(threshBinary,3)
    shapeProps = regionprops("table", threshBinary(:,:,fr),"Area", "BoundingBox","Centroid", "SubarrayIdx","PixelIdxList");
    shapeProps = sortrows(shapeProps,"Area", "descend");

    % limit to objects larger than min
    rows = shapeProps.Area > waveMinSize;

    if sum(rows) > 0
        shapeProps = shapeProps(rows,:);
        for a = 1:height(shapeProps)
            shapeProps.Frame(a) = fr;
        end
        waveTable = [waveTable; shapeProps];
    end
end

% clean for duplicates
[~,indx] = unique(waveTable.Centroid(:,1),'rows');
waveTable = waveTable(indx,:);
waveTable = sortrows(waveTable,"Frame", "ascend");


% group by bounding box overlap
frames = unique(waveTable.Frame);
for fr = frames'
    frameObjIn = waveTable.Frame == fr;
    currentFrame = waveTable(frameObjIn,:);

    if height(currentFrame) > 1
        currentFrameBB = waveTable.BoundingBox(frameObjIn,:);
        currentAreas = waveTable.Area(frameObjIn,:);

        % get possible combinations for overlap
        possCombs = nchoosek(1:size(currentFrameBB,1), 2);
        possCombs(:,3) = 0;

        % check overlap
        for cmbs = 1:size(possCombs,1)
            overlap = bboxOverlapRatio(currentFrameBB(possCombs(cmbs,1),:), currentFrameBB(possCombs(cmbs,2),:));

            if overlap > 0
                % find larger area for bounding box
                [~, indMax] = max([currentAreas(possCombs(cmbs,1)) currentAreas(possCombs(cmbs,2))]);

                groupBB = possCombs(cmbs,indMax);
                possCombs(cmbs,3) = groupBB;
            end
        end

        % remove non overlapping combinations
        indices = find(possCombs(:,3)==0);
        possCombs(indices,:) =[];

        % get the indices of BBs
        indexBB = find(frameObjIn);

        % copies the grouping BB into the waveTable;
        for ss = 1:size(possCombs,1)
            waveTable.Group(indexBB(possCombs(ss,1))) = possCombs(ss,3);
            waveTable.Group(indexBB(possCombs(ss,2))) = possCombs(ss,3);
        end

%         imshow(tifTreshFilled(:,:,fr))
%         hold on
%         for vv = 1:size(currentFrameBB,1)
%             rectangle('Position', currentFrameBB(vv,:), 'EdgeColor','r');
%         end
%         title(['Frame No: ' num2str(fr)]);

    else
        waveTable.Group(frameObjIn) = 0;
    end
end

disp('Finished getting Ca objects')

%% check for bounding box overlap across frames and start wave categorization

% add overlap BB index coloumn
waveTable.OverlapBBindx = cell(height(waveTable),1);

% for each frame
for rr = 1:size(threshBinary,3)-1
    % check if there are wave objects in that frame
    if ismember(rr, waveTable.Frame(:))
        % get frame objects
        frameObjIn = waveTable.Frame == rr;
        frameIndNum = find(frameObjIn);
        currentFrame = waveTable(frameObjIn,:);

        % get next frame objects
        frameObjInNext = waveTable.Frame == rr+1;

        % if there are any objects to match in the next frame
        if sum(frameObjInNext) > 0
            nextFrame = waveTable(frameObjInNext,:);


            % for current frame objects
            counter = 0;
            for currentBB = 1:height(currentFrame)

                % for all next frame objects mark the overlap with the
                % previous frame object indexes
                for nextFrameBB = 1:height(nextFrame)
                    overlap = bboxOverlapRatio(currentFrame.BoundingBox(currentBB,:), nextFrame.BoundingBox(nextFrameBB,:));

                    if overlap > 0
                        counter = counter +1;
                        nextFrame.OverlapBBindx{nextFrameBB}(counter) = frameIndNum(currentBB);
                    end
                end
            end

            % find and remove all zeros (small bug)
            idxZeros = cellfun(@(c)(ismember(c,0)), nextFrame.OverlapBBindx, 'UniformOutput',false);

            for u = 1:height(nextFrame)
                if ~isempty(nextFrame.OverlapBBindx{u}) && sum(nextFrame.OverlapBBindx{u}) > 0
                    nextFrame.OverlapBBindx{u} = nextFrame.OverlapBBindx{u}(~idxZeros{u});
                end
            end

            % overwrite frame object info
            waveTable(frameObjInNext,:) = nextFrame;
        else

        end
    end
end

disp('Checking for overlaps');
%% classify objects into waves
waveTable.waveNumber = zeros(height(waveTable),1);
currentWave = 0;

for ob = 1:height(waveTable)

    if isempty(waveTable.OverlapBBindx{ob})
        currentWave = currentWave +1;
        waveTable.waveNumber(ob) = currentWave;
    else
        % get the minimum wave number
        minWaveNo = min(waveTable.waveNumber(waveTable.OverlapBBindx{ob}));

        % overwrite the grouped previous objects
        waveTable.waveNumber(waveTable.OverlapBBindx{ob}) = minWaveNo;
        waveTable.waveNumber(ob) = minWaveNo;

    end
end

% renumber the wave nums
[~,~,newWaveNo] = unique(waveTable.waveNumber);
waveTable.waveNumber = newWaveNo;


%% add wave table to waves struct
waves.waveTable = waveTable;

%% colorize imagestack to check waves
waveCols = distinguishable_colors(length(unique(waveTable.waveNumber)),{'w','k'});

% reshape into RGB image stack
tifGauSubRGB = repmat(tifTreshFilled, 1,1,1,3);
frames = unique(waveTable.Frame);

% for frames with waves
for fr = frames'

    % get frame objects
    frameObjIn = waveTable.Frame == fr;
    currentFrame = waveTable(frameObjIn,:);

    % for all objects in frame, colorize pixels
    for BB = 1:height(currentFrame)
        currentFrameX = currentFrame.SubarrayIdx{BB,1} ;
        currentFrameY = currentFrame.SubarrayIdx{BB,2} ;

        currentFramePix = currentFrame.PixelIdxList{BB};

        % initialise RGB
        R = tifGauSubRGB(:,:,fr,1);
        G = tifGauSubRGB(:,:,fr,2);
        B = tifGauSubRGB(:,:,fr,3);

        R(currentFramePix) = R(currentFramePix) * waveCols(currentFrame.waveNumber(BB),1);
        G(currentFramePix) = G(currentFramePix) * waveCols(currentFrame.waveNumber(BB),2);
        B(currentFramePix) = B(currentFramePix) * waveCols(currentFrame.waveNumber(BB),3);

        tifGauSubRGB(:,:,fr,:) = cat(3, R,G,B);
    end
end

% save RGB timeseries stack
options.color = 1;
saveastiff(permute(tifGauSubRGB, [1 2 4 3]), fullfile(folderPath, [name(1:end-5) '_waveCol2.tif']), options);

%% create wave extent images

% get image to color over
SDImage = std(double(tifStack),[],3);
SDImage = imadjust(rescale(SDImage));

% color over each wave extent
for w = 1:max(waveTable.waveNumber)

    % get all objects for a single waves
    subTabIndx = waveTable.waveNumber == w;
    subTable = waveTable(subTabIndx,:);

    % get all the pixels involved in wave
    wavePixels = unique(cat(1,subTable.PixelIdxList{:}));
    [wavePixX, wavePixY] = ind2sub(size(SDImage), wavePixels);

    % get wave extent in pixels
    waves.waveArea(w) = length(wavePixels);
    waves.waveAreaMicron(w) = length(wavePixels) * exStruct.downsampledRes;


    SDImageRGB = repmat(SDImage, 1, 1, 3);

    % colorize each pixel
    for pix = 1:length(wavePixX)
    SDImageRGB(wavePixX(pix), wavePixY(pix) ,:) = squeeze(SDImageRGB( wavePixX(pix), wavePixY(pix),:))' .* waveCols(w,:);
    end

    % save the things
    wavePicDir = fullfile(folderPath, [name(1:end-5) '_waveMaps']);

    if ~exist(wavePicDir)
        mkdir(wavePicDir);
    end

    imwrite(SDImageRGB, fullfile(wavePicDir, sprintf('%s_wave_%03d.tif', name(1:end-5), w)));
end

%% get the wave metrics

for w = 1:max(waveTable.waveNumber)

    indxWave = waveTable.waveNumber ==w;
    currentWave = waveTable(indxWave,:);

    % wave frames
    waveFrameFirst(w) = currentWave.Frame(1);
    waveFrameLast(w) = currentWave.Frame(end);

    % wave time on/off
    waveTimeOn(w) = exStruct.frameInfo.frameTime(waveFrameFirst(w));
    waveTimeOff(w) = exStruct.frameInfo.frameTime(waveFrameLast(w));

    % average DF/F
    DF_PixelList = [];
    for cc = 1:height(currentWave)
        currentDF_Frame = exStruct.dF(:,:,currentWave.Frame(cc));
        currentPixels = currentDF_Frame(currentWave.PixelIdxList{cc});
        DF_PixelList = [DF_PixelList; currentPixels ];
    end

    waveDFAverage(w) = mean(DF_PixelList);
end

%% add into waves
waves.waveFrameFirst = waveFrameFirst;
waves.waveFrameLast = waveFrameLast;
waves.waveTimeOn = waveTimeOn;
waves.waveTimeOff = waveTimeOff;
waves.waveDFAverage = waveDFAverage;

exStruct.waves = waves;

save(fullfile(folderPath,[name(1:end-5) '_ExStruct.mat']), 'exStruct');

end