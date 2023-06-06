function checkWaveOrigins(exStructPath)

%% load in exStruct
exStruct = load(exStructPath);
exStruct = exStruct.exStruct;

waves = exStruct.waves;

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
