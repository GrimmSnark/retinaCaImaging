function trackCalciumWaves(filepathDF)

%% defaults
waveMinSize = 300;

%% load in dF movie
tifStack = read_Tiffs(filepathDF);
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


%% threshold
tifThresh = imbinarize(tifGauSub, 0.01);
tifThresh = im2uint8(tifThresh);


%% fill holes to make wave detection easier
for cc = 1:size(tifThresh,3)
    tifTreshFilled(:,:,cc) = imfill(tifThresh(:,:,cc));
end

threshBinary = imbinarize(tifTreshFilled);

% %% optic flow
% opticFlow =opticalFlowFarneback;
% for ee = 1:size(tifGauSub,3) -1
%     flow = estimateFlow(opticFlow, double(threshBinary(:,:,ee)));
%     imshow(imadjust(tifGauSub(:,:,ee)));
%     hold on
%     plot(flow,'DecimationFactor',[5 5],'ScaleFactor',1);
%     hold off
%     pause
%
% end
%

%%
waveTable = [];
opticFlow =opticalFlowHS;
for fr = 1:size(threshBinary,3)
    shapeProps = regionprops("table", threshBinary(:,:,fr),"Area", "BoundingBox","Centroid", "SubarrayIdx","PixelIdxList");
    shapeProps = sortrows(shapeProps,"Area", "descend");

    % limit to objects larger than min
    rows = shapeProps.Area > waveMinSize;

    if sum(rows) > 0
        shapeProps = shapeProps(rows,:);
        for a = 1:height(shapeProps)
            flow = estimateFlow(opticFlow, double(tifGauSub(:,:,fr)));
            shapeProps.Frame(a) = fr;

            % get the optic flow stuff
            pixelIndxs = shapeProps.PixelIdxList{a};

            % angles
            pixelAngles = rad2deg(flow.Orientation(pixelIndxs));
            % correct for negative angles
            pixelAngles(pixelAngles< 0) = pixelAngles(pixelAngles< 0) + 360;
            % magnitudes
            pixelMag = flow.Magnitude(pixelIndxs);

            % vector average
            meanVec = mean_vector_direction_magnitude([pixelAngles pixelMag]);

            % add into shapeProps
            shapeProps.meanOrientation(a) = meanVec.mean_angle_degrees;
            shapeProps.meanMagnitude(a) = mean(pixelMag);
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

    else
        waveTable.Group(frameObjIn) = 0;
    end
end

%% check for bounding box overlap across frames and start wave categorization

% add overlap BB index coloumn
waveTable.OverlapBBindx = zeros(height(waveTable),1);

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
            for currentBB = 1:height(currentFrame)

                % for all next frame objects
                for nextFrameBB = 1:height(nextFrame)
                    overlap = bboxOverlapRatio(currentFrame.BoundingBox(currentBB,:), nextFrame.BoundingBox(nextFrameBB,:));

                    if overlap > 0
                        nextFrame.OverlapBBindx(nextFrameBB) = frameIndNum(currentBB);
                    end
                end
            end

            waveTable(frameObjInNext,:) = nextFrame;
        else

        end
    end
end

%% classify objects into waves
waveTable.waveNumber = zeros(height(waveTable),1);
currentWave = 0;

for ob = 1:height(waveTable)

    if waveTable.OverlapBBindx(ob) == 0
        currentWave = currentWave +1;
        waveTable.waveNumber(ob) = currentWave;
    elseif waveTable.OverlapBBindx(ob) > 0
        waveTable.waveNumber(ob) = waveTable.waveNumber(waveTable.OverlapBBindx(ob));
    end
end

%% colorize imagestack to check waves
waveCols = distinguishable_colors(length(unique(waveTable.waveNumber)),{'w','k'});

% reshape into RGB image stack
tifGauSubRGB = permute(repmat(tifTreshFilled, 1,1,1,3), [1 ,2 ,4 ,3]);
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
        R = tifGauSubRGB(:,:,1,fr);
        G = tifGauSubRGB(:,:,2,fr);
        B = tifGauSubRGB(:,:,3,fr);

        R(currentFramePix) = R(currentFramePix) * waveCols(currentFrame.waveNumber(BB),1);
        G(currentFramePix) = G(currentFramePix) * waveCols(currentFrame.waveNumber(BB),2);
        B(currentFramePix) = B(currentFramePix) * waveCols(currentFrame.waveNumber(BB),3);

        tifGauSubRGB(:,:,:,fr) = cat(3, R,G,B);
    end
end

[folderPath, name] = fileparts(filepathDF);
% options.color = 1;
% saveastiff(tifGauSubRGB, fullfile(folderPath, [name(1:end-5) '_waveCol.tif']), options);

% save RGB timeseries stack
bfsave(tifGauSubRGB, fullfile(folderPath, [name(1:end-5) '_waveCol.tif']), 'dimensionOrder', 'XYCTZ' );

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
    waveArea(w) = length(wavePixels);

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

%% get average dF/F for each waves

for w = 1:max(waveTable.waveNumber)

end
end