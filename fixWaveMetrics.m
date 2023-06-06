function fixWaveMetrics(exStructPath)


exStruct = load(exStructPath);
exStruct = exStruct.exStruct;

waveTable = exStruct.waves.waveTable;

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

%% create wave extent images

% color over each wave extent
for w = 1:max(waveTable.waveNumber)

    % get all objects for a single waves
    subTabIndx = waveTable.waveNumber == w;
    subTable = waveTable(subTabIndx,:);

    % get all the pixels involved in wave
    wavePixels = unique(cat(1,subTable.PixelIdxList{:}));

    % get wave extent in pixels
    waves.waveArea(w) = length(wavePixels);
    waves.waveAreaMicron(w) = length(wavePixels) * exStruct.downsampledRes;
end

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


exStruct.waves = waves;

%%
%% run through waves
waveStartInCl = false(size(waves.waveArea))';
wavePassInCl = false(size(waves.waveArea))';
edgeDistanceWaveStartPix = nan(size(waves.waveArea))';
edgeDistanceWaveStartPix = nan(size(waves.waveArea))';


for i = 1:length(waves.waveArea)
    waveObjectsIndx = waves.waveTable.waveNumber == i;
    waveObjects = waves.waveTable(waveObjectsIndx,:);

    [~, minFrameInd] = min(waveObjects.Frame);

    minObjectHull = waveObjects.ConvexHull(minFrameInd);

    for cl = 1:height(exStruct.clusterCells.clusterPropsDS)
        % check if the wave starts in the cluster
        for ob = 1:length(minObjectHull)
            clusterHull = exStruct.clusterCells.clusterPropsDS.ConvexHull{cl};
            in = polyxpoly(minObjectHull{ob}(:,1), minObjectHull{ob}(:,2),clusterHull(:,1), clusterHull(:,2));

            minObjectPoly.x = minObjectHull{ob}(:,1);
            minObjectPoly.y = minObjectHull{ob}(:,2);

            clusterP.x = clusterHull(:,1);
            clusterP.y =  clusterHull(:,2);
            edgeDistanceWaveStartPix(i,cl) = min_dist_between_two_polygons(minObjectPoly, clusterP);

            % plots the wave objects and cluster locations

            %             mask = exStruct.clusterCells.maskDS;
            %             mask(waveObjects.PixelIdxList{minFrameInd}) = 1;
            %             imshow(mask);
            %             hold on
            %             mapshow(minObjectHull{ob}(:,1), minObjectHull{ob}(:,2), 'Marker','+');
            %             mapshow(clusterHull(:,1), clusterHull(:,2), 'Marker','+', 'Color','r');
            %             title(['Cluster No: ' num2str(cl) ' Wave No: ' num2str(i) ' Frame No: ' num2str(waveObjects.Frame(minFrameInd))]);
            %             xlim([0 512]);
            %             ylim([0 512]);
            %             close(gcf)

            if ~isempty(in)
                waveStartInCl(i) = true;
            end

            % check if any wave object passes through cluster
            for waveOb = 1:height(waveObjects)
                in = polyxpoly(waveObjects.ConvexHull{waveOb}(:,1), waveObjects.ConvexHull{waveOb}(:,2),clusterHull(:,1), clusterHull(:,2));

                %                 mapshow(waveObjects.ConvexHull{waveOb}(:,1), waveObjects.ConvexHull{waveOb}(:,2), 'Marker','+');
                %                 hold on
                %                 mapshow(clusterHull(:,1), clusterHull(:,2), 'Marker','+', 'Color','r');
                %                 title(['Cluster No: ' num2str(cl) ' Wave No: ' num2str(i) ' Frame No: ' num2str(waveObjects.Frame(waveOb))]);
                %                 xlim([0 2048]);
                %                 ylim([0 2048]);
                %                 close(gcf)


                if ~isempty(in)
                    wavePassInCl(i) = true;
                end
            end
        end
    end
end

%% save stuff into exStruct
exStruct.waves.waveStartInCl = waveStartInCl;
exStruct.waves.wavePassInCl = wavePassInCl;
exStruct.waves.edgeDistanceWaveStartClPix = edgeDistanceWaveStartPix;
exStruct.waves.edgeDistanceWaveStartClMicron = edgeDistanceWaveStartPix/exStruct.image.pixelSize;





% save data
save(exStructPath, "exStruct", '-v7.3');


end