function [waves,colStack] = recalculateWaveMetrics(app)

waveTable = app.exStruct.waves.waveTable;
exStruct = app.exStruct;
waves = app.exStruct.waves;

%% get the wave metrics
waveDFAverage = [];
for w = 1:max(waveTable.waveNumber)

    indxWave = waveTable.waveNumber ==w;
    currentWave = waveTable(indxWave,:);

    % wave frames
    waveFrameFirst(w) = currentWave.Frame(1);
    waveFrameLast(w) = currentWave.Frame(end);

    % wave time on/off
    waveTimeOn(w) = exStruct.frameInfo.frameTime(waveFrameFirst(w));
    waveTimeOff(w) = exStruct.frameInfo.frameTime(waveFrameLast(w));

    % wave DF average (not really useful)
    try
        DF_PixelList = [];
        for cc = 1:height(currentWave)
            currentDF_Frame = exStruct.dF(:,:,currentWave.Frame(cc));
            currentPixels = currentDF_Frame(currentWave.PixelIdxList{cc});
            DF_PixelList = [DF_PixelList; currentPixels ];
        end

        waveDFAverage(w) = mean(DF_PixelList);
    catch

    end

    % wave trajectory
    count = 1;
    for fr = waveFrameLast(w):-1: waveFrameFirst(w)
        currentWaveFrame = currentWave(currentWave.Frame == fr,:);

        if ~isempty(currentWaveFrame)

            % if more than one frame object find closest
            if height(currentWaveFrame) > 1 && count > 1 % if in the middle of the wave

                % weighted distance by size of object
                [objectDis  ] = pdist2(currentWaveFrame.Centroid, centerPerFrame{w}(count-1,:), 'euclidean');
                weightedDistance = currentWaveFrame.Area ./ objectDis;
                [~, in] = max(weightedDistance);

            elseif height(currentWaveFrame) > 1 && count == 1 % if at the start of the wave with mutiple objects

                [~, in] = max(currentWaveFrame.Area);
            else

                in = 1;
            end

            centerPerFrame{w}(count,:) = currentWaveFrame.Centroid(in,:);
            count = count + 1;
        end
    end

    % flip the trajectory to the start
    centerPerFrame{w} = flip( centerPerFrame{w});

    if length(centerPerFrame{w}) >2
        for i = 1:length(centerPerFrame{w})-1
            % centroid distance per frame
            distancePerFramePix{w}(i)= pdist([centerPerFrame{w}(i+1,:) ;centerPerFrame{w}(i,:)], 'euclidean');
            distancePerFrameMicron{w}(i)= distancePerFramePix{w}(i) * exStruct.image.pixelSize;

            % speed per frame
            speedPerFrameMicronSec{w}(i) = distancePerFrameMicron{w}(i) / exStruct.framePeriod;
            maxSpeed(w) = max(speedPerFrameMicronSec{w});
        end
    end
end

%% colorize imagestack to check waves
waveCols = distinguishable_colors(length(unique(waveTable.waveNumber)),{'w','k'});

% rescale to 0-1
colStack = rescale(app.disStack.waveColStack);

% find max across RGB
colStackMax = max(colStack,[], 4);
colStack = repmat(colStackMax,1,1,1,3); % recreate colStack so that any channels which are black, ie RGB get turned to white in line below

colStack(colStack > 0) = 1;

% reshape into RGB image stack
frames = unique(waveTable.Frame);

% for frames with waves
for fr = frames'

    % get frame objects
    frameObjIn = waveTable.Frame == fr;
    currentFrame = waveTable(frameObjIn,:);

    % for all objects in frame, colorize pixels
    for BB = 1:height(currentFrame)
%         currentFrameX = currentFrame.SubarrayIdx{BB,1} ;
%         currentFrameY = currentFrame.SubarrayIdx{BB,2} ;

        currentFramePix = currentFrame.PixelIdxList{BB};

        % initialise RGB
        R = colStack(:,:,fr,1);
        G = colStack(:,:,fr,2);
        B = colStack(:,:,fr,3);

        R(currentFramePix) = R(currentFramePix) * waveCols(currentFrame.waveNumber(BB),1);
        G(currentFramePix) = G(currentFramePix) * waveCols(currentFrame.waveNumber(BB),2);
        B(currentFramePix) = B(currentFramePix) * waveCols(currentFrame.waveNumber(BB),3);

        colStack(:,:,fr,:) = cat(3, R,G,B);
    end
end

% save RGB timeseries stack
% options.color = 1;
% colStack = im2uint8(colStack);
% saveastiff(permute(colStack, [1 2 4 3]), fullfile([exStruct.filePath(1:end-4) '_waveCol.tif']), options);

%% create wave extent images

% % get image to color over
% SDImage = std(double(app.disStack.dFTifStack),[],3);
% SDImage = imadjust(rescale(SDImage));
% 
% % color over each wave extent
% for w = 1:max(waveTable.waveNumber)
% 
%     % get all objects for a single waves
%     subTabIndx = waveTable.waveNumber == w;
%     subTable = waveTable(subTabIndx,:);
% 
%     % get all the pixels involved in wave
%     wavePixels = unique(cat(1,subTable.PixelIdxList{:}));
%     [wavePixX, wavePixY] = ind2sub(size(SDImage), wavePixels);
% 
%     % get wave extent in pixels
%     waves.waveArea(w) = length(wavePixels);
%     waves.waveAreaMicron(w) = length(wavePixels) * exStruct.downsampledRes;
% 
% 
%     SDImageRGB = repmat(SDImage, 1, 1, 3);
% 
%     % colorize each pixel
%     for pix = 1:length(wavePixX)
%         SDImageRGB(wavePixX(pix), wavePixY(pix) ,:) = squeeze(SDImageRGB( wavePixX(pix), wavePixY(pix),:))' .* waveCols(w,:);
%     end
% 
%     % save the things
%     wavePicDir = [exStruct.filePath(1:end-4) '_waveMaps'];
% 
%     if ~exist(wavePicDir)
%         mkdir(wavePicDir);
%     else
%         if w == 1
%             rmdir(wavePicDir, 's');
%             mkdir(wavePicDir);
%         end
%     end
% 
%     [~, fileName] = fileparts(exStruct.filePath);
% 
%     imwrite(SDImageRGB, fullfile(wavePicDir, sprintf('%s_wave_%03d.tif', fileName, w)));
% end
% 

%% add into waves
waves.waveTable = waveTable;
waves.waveFrameFirst = waveFrameFirst;
waves.waveFrameLast = waveFrameLast;
waves.waveTimeOn = waveTimeOn;
waves.waveTimeOff = waveTimeOff;
waves.waveDFAverage = waveDFAverage;
waves.centerPerFrame = centerPerFrame;
waves.distancePerFramePix = distancePerFramePix;
waves.distancePerFrameMicron = distancePerFrameMicron;
waves.speedPerFrameMicronSec = speedPerFrameMicronSec;


end