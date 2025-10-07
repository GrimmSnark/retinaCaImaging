function calculateYoPro_D1_2(folderPath)
% Calculates the D1/D2 ratio for yopro-1 stained retinal
% wholemounts. This will run on a full folder of images. Requires each
% image in the folder to be processed for to a binary image with cells in
% black

%% open FIJI
% initalize MIJI
intializeMIJ;

RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();
RC.reset();

%% find the images we want to use ie 'mask'
% maskFilepath = dir([folderPath '\*mask2*']);
%
% if isempty(maskFilepath)
%     maskFilepath = dir([folderPath '\*Mask2*']);
% end

maskFilepath = dir([folderPath '**\*mask*']);

if isempty(maskFilepath)
    maskFilepath = dir([folderPath '**\*Mask*']);
end


for x = 1:length(maskFilepath)
        filePathTemp = maskFilepath(x);

        pattern = '[pP](\d+)';
        pDay = regexp(filePathTemp.name, pattern, 'tokens');
        pDay = str2double(pDay{1}{:});

        filePathTemp = fullfile(filePathTemp.folder, filePathTemp.name);

        masks = read_Tiffs(filePathTemp);
        % masks = ~masks;
        yoproCells = bwconncomp(masks);
        yoproLabelled = labelmatrix(yoproCells);
        pixelImSize = yoproCells.ImageSize;
        centroidStruct = regionprops(yoproLabelled,"Centroid");
        maskCentroid = vertcat(centroidStruct.Centroid);

        %% Get the optic nerve head and retinal bounds

        if exist([filePathTemp(1:end-9) '_ROIs.zip'])
            RC.open([filePathTemp(1:end-9) '_ROIs.zip']);
        elseif exist([filePathTemp(1:end-9) 'ROIs.zip'])
            RC.open([filePathTemp(1:end-9) 'ROIs.zip']);
        elseif exist([filePathTemp(1:end-18) '_ROIs.zip'])
            RC.open([filePathTemp(1:end-18) '_ROIs.zip']);
        end

        % RC.open([filePathTemp(1:end-18) '_ROIs.zip']);

        ROIobjects = RC.getRoisAsArray;

        % optic nerve poly
        opticNerveMask = createLabeledROIFromImageJPixels(yoproCells.ImageSize ,ROIobjects(1));
        centroidStruct = regionprops(opticNerveMask,"Centroid");
        opticNerveCentroid = centroidStruct.Centroid;

        % retina boundary poly
        retinaBoundMask = createLabeledROIFromImageJPixels(yoproCells.ImageSize  ,ROIobjects(2));
        retinaBoundPoly = bwboundaries(retinaBoundMask');
        retinaBoundShape = polyshape(retinaBoundPoly{1});


        %%

        distFromCenter = pdist2(opticNerveCentroid,maskCentroid);
        count = 1;
        center2RetinaEdge = [];
        retinaEdgePos = [];
        pFit = [];
        for w = 1:length(maskCentroid)
            pFit(w,:,:) = fitStraightLine(opticNerveCentroid, maskCentroid(w,:), [0 max(pixelImSize)]);

            % get the points in and out of the retina shape
            currLine = squeeze(pFit(w,:,:));
            [inR, outR] = intersect(retinaBoundShape,  currLine);

            % plot(retinaBoundShape)
            % hold on
            % scatter(maskCentroid(w,1), maskCentroid(w,2));
            % scatter(opticNerveCentroid(1), opticNerveCentroid(2));
            % plot(currLine(:,1), currLine(:,2));
            % close;

            % error catch for objects outside of retina bounds
            if isempty(inR)
                continue
            end
            % %%%%

            [retBoundDist,indUsed] = pdist2(inR, maskCentroid(w,:),'euclidean','Smallest',1);
            line2Use = [opticNerveCentroid; inR(indUsed,:) ];

            center2RetinaEdge(count) = pdist2(inR(indUsed,:), opticNerveCentroid);
            retinaEdgePos(count,:) = inR(indUsed,:);

            count = count +1;
            % use this to check calculations

            % g = imshow(imbinarize(masks));
            % hold on
            % plot(retinaBoundShape)
            % scatter(opticNerveCentroid(1), opticNerveCentroid(2));
            % scatter(maskCentroid(w,1), maskCentroid(w,2));
            %
            % plot(line2Use(:,1),line2Use(:,2),'b');

        end

        distFromCenter = distFromCenter';
        center2RetinaEdge =center2RetinaEdge';

        %% build table
        relDistanceTab = table;
        relDistanceTab = table(maskCentroid, distFromCenter, pFit, center2RetinaEdge , retinaEdgePos);
        relDistanceTab.pDay = repmat(pDay,height(relDistanceTab),1);
        relDistanceTab.D1_2 = relDistanceTab.distFromCenter./relDistanceTab.center2RetinaEdge;

        % clean numbers
        relDistanceTab(relDistanceTab.D1_2>1,:) = [];
        relDistanceTab(relDistanceTab.D1_2<0,:) = [];

        D1_2Table = table();
        D1_2Table.pDay = relDistanceTab.pDay;
        D1_2Table.D1_2 = relDistanceTab.D1_2;

        %% save as excel

        writetable(D1_2Table, [filePathTemp(1:end-9) '_D1_2_v2.xlsx']);

    RC.close
end
end