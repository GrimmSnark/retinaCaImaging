function calculateHMXO_D1_2(folderPath)
% Calculates the D1/D2 ratio for HMOX-1 microglia stained retinal
% wholemounts. Requires each image in the folder to be processed for to a
% binary image with cells in black


%% open FIJI
% initalize MIJI
intializeMIJ;

RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();
RC.reset();

%% find the images we want to use ie 'mask'
maskFilepath = dir([folderPath '\*mask*']);

if isempty(maskFilepath)
    maskFilepath = dir([folderPath '\*Mask*']);
end

maskFilepath = maskFilepath(1);

pattern = '[pP](\d+)';
pDay = regexp(maskFilepath.name, pattern, 'tokens');
pDay = str2double(pDay{1}{:});

maskFilepath = fullfile(maskFilepath.folder, maskFilepath.name);

masks = read_Tiffs(maskFilepath);
masks = imbinarize(masks);
pixelIm = size(masks,1);
centroidStruct = regionprops(masks,"Centroid"); 
maskCentroid = vertcat(centroidStruct.Centroid);

%% Get the optic nerve head and retinal bounds
roiPath = dir([folderPath '\*.zip']);
roiPath = fullfile(roiPath.folder, roiPath.name);

RC.open(roiPath);
ROIobjects = RC.getRoisAsArray;

% optic nerve poly
opticNerveMask = createLabeledROIFromImageJPixels([pixelIm pixelIm] ,ROIobjects(1));
centroidStruct = regionprops(opticNerveMask,"Centroid");
opticNerveCentroid = centroidStruct.Centroid;

% retina boundary poly
retinaBoundMask = createLabeledROIFromImageJPixels([pixelIm pixelIm] ,ROIobjects(2));
retinaBoundPoly = bwboundaries(retinaBoundMask');
retinaBoundShape = polyshape(retinaBoundPoly{1});


%%

distFromCenter = pdist2(opticNerveCentroid,maskCentroid);

for w = 1:length(maskCentroid)
    pFit(w,:,:) = fitStraightLine(opticNerveCentroid, maskCentroid(w,:), [0 pixelIm]);

    % get the points in and out of the retina shape
    currLine = squeeze(pFit(w,:,:));
    [inR, outR] = intersect(retinaBoundShape,  currLine);
    [retBoundDist,indUsed] = pdist2(inR, maskCentroid(w,:),'euclidean','Smallest',1);
    line2Use = [opticNerveCentroid; inR(indUsed,:) ];

    center2RetinaEdge(w) = pdist2(inR(indUsed,:), opticNerveCentroid);
    retinaEdgePos(w,:) = inR(indUsed,:);

    % use this to check calculations
%     g = imshow(imadjust(masks));
%     hold on
%     plot(retinaBoundShape)
%     scatter(opticNerveCentroid(1), opticNerveCentroid(2));
%     scatter(maskCentroid(w,1), maskCentroid(w,2));
%     plot(line2Use(:,1),line2Use(:,2),'b');

end

distFromCenter = distFromCenter';
center2RetinaEdge =center2RetinaEdge';

%% build table
relDistanceTab = table;
relDistanceTab = table(maskCentroid, distFromCenter, pFit, center2RetinaEdge , retinaEdgePos);
relDistanceTab.pDay = repmat(pDay,height(relDistanceTab),1);
relDistanceTab.D1_2 = relDistanceTab.distFromCenter./relDistanceTab.center2RetinaEdge;

D1_2Table = table();
D1_2Table.pDay = relDistanceTab.pDay;
D1_2Table.D1_2 = relDistanceTab.D1_2;

%% save as excel

writetable(D1_2Table, [folderPath '\D1_2.xlsx']);
end