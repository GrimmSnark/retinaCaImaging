function plotRelativeDisWaveStarts(exStructPath)

%% load in exStruct
exStruct = load(exStructPath);
disp('Loaded in exStruct.mat')

try
    exStruct = exStruct.exStruct;
catch
    exStruct = exStruct.metaData;
end

%% open FIJI
% initalize MIJI
intializeMIJ;

RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();
RC.reset();

%% load in images
SD_imagePath = dir([exStructPath(1:end-13) '_dF_SD.tif']);
SD_image = read_Tiffs(fullfile(SD_imagePath.folder,SD_imagePath.name));


if ~isfield(exStruct.waves, 'relDistanceTab')
    %% load into FIJI
    stackImpSD = MIJ.createImage('SD', SD_image, 1);

    ij.IJ.runMacro('setTool("polygon");');

    %% choose ROI shapes for cluster cells

    % Sets up diolg box to allow for user input to choose cell ROIs
    opts.Default = 'Continue';
    opts.Interpreter = 'tex';

    questText = [{'Draw ROIs around the optic nerve head and retina'} ...
        {'Add ROI to ROI Manager by "t" for optic nerve head then retina'} {''} ...
        {'If you are happy to move on with analysis click  \bfContinue\rm'} ...
        {'Or click  \bfExit Script\rm or X out of this window to exit script'}];

    response = questdlg(questText, ...
        'Check and choose ROIs', ...
        'Continue', ...
        'Exit Script', ...
        opts);


    % deals with repsonse
    switch response

        case 'Continue' % if continue, goes on with analysis
            ROInumber = RC.getCount();
            disp(['You have selected ' num2str(ROInumber) ' ROIs, moving on...']);

        case 'Exit Script' % if you want to exit and end
            killFlag = 1;
            return
        case ''
            killFlag = 1; % if you want to exit and end
            return
    end

    %% create ROI image
    ROIobjects = RC.getRoisAsArray;
    pixelIm = stackImpSD.getWidth;

    % optic nerve poly
    opticNerveMask = createLabeledROIFromImageJPixels([pixelIm pixelIm] ,ROIobjects(1));
    centroidStruct = regionprops(opticNerveMask,"Centroid");
    opticNerveCentroid = centroidStruct.Centroid;

    % retina boundary poly
    retinaBoundMask = createLabeledROIFromImageJPixels([pixelIm pixelIm] ,ROIobjects(2));
    retinaBoundPoly = bwboundaries(retinaBoundMask');
    retinaBoundShape = polyshape(retinaBoundPoly{1});

else
    % optic nerve poly
    centroidStruct = regionprops(exStruct.waves.opticNerveMask,"Centroid");
    opticNerveCentroid = centroidStruct.Centroid;

    % retina boundary poly
    retinaBoundPoly = bwboundaries(exStruct.waves.retinaBoundMask');
    retinaBoundShape = polyshape(retinaBoundPoly{1});
end


%% get distances from optic nerve head, retina boundary, and blood plexus

for x = 1:length(exStruct.waves.centerPerFrame)
    waveStarts(x,:) = exStruct.waves.centerPerFrame{x}(1,:);
end

distFromCenter = pdist2(opticNerveCentroid,waveStarts);


BV_shape = polyshape(exStruct.waves.BV_poly);
% fit line from center to wave start point and calcuate distances
for w = 1:length(exStruct.waves.centerPerFrame)
    pFit(w,:,:) = fitStraightLine(opticNerveCentroid, waveStarts(w,:), [0 512]);

    % get the points in and out of the retina shape
    currLine = squeeze(pFit(w,:,:));
    [inR, outR] = intersect(retinaBoundShape,  currLine);
    [retBoundDist,indUsed] = pdist2(inR, waveStarts(w,:),'euclidean','Smallest',1);
    line2Use = [opticNerveCentroid; inR(indUsed,:) ];

    center2RetinaEdge(w) = pdist2(inR(indUsed,:), opticNerveCentroid);
    retinaEdgePos(w,:) = inR(indUsed,:);

    % get the points in and out of the BV
    [inBV, outBV] = intersect(BV_shape,  line2Use);

    % use point not equal to the retina center position
    index = and(any(inBV == opticNerveCentroid, 2), any(inBV == opticNerveCentroid, 2));

    % use point not equal to the retina center position
    index2Del = and(any(inBV == opticNerveCentroid, 2), any(inBV == opticNerveCentroid, 2));
    inBV(index2Del,:) = [];


    [BV_BoundDist,indUsed] = pdist2(inBV, waveStarts(w,:),'euclidean','Smallest',1);
    line2UseBV = [opticNerveCentroid; inBV(indUsed,:) ];

    center2BVEdge(w) = pdist2(inBV(indUsed,:), opticNerveCentroid);
    BVEdgePos(w,:) = inBV(indUsed,:);

    % use this to check calculations
%     g = imshow(imadjust(SD_image));
%     hold on
%     plot(retinaBoundShape)
%     plot(BV_shape)
%     scatter(opticNerveCentroid(1), opticNerveCentroid(2));
%     scatter(waveStarts(w,1), waveStarts(w,2));
%     plot(line2Use(:,1),line2Use(:,2),'b',line2UseBV(:,1),line2UseBV(:,2),'r');

end

distFromCenter = distFromCenter';
center2RetinaEdge =center2RetinaEdge';
center2BVEdge = center2BVEdge';

%% build table
relDistanceTab = table(waveStarts, distFromCenter, pFit, center2RetinaEdge , retinaEdgePos, center2BVEdge, BVEdgePos);
relDistanceTab.D1_2 = relDistanceTab.distFromCenter./relDistanceTab.center2RetinaEdge;
relDistanceTab.D1_3 = relDistanceTab.distFromCenter./relDistanceTab.center2BVEdge;

%% add to exStruct and save

if ~isfield(exStruct.waves, 'relDistanceTab')
    exStruct.waves.relDistanceTab = relDistanceTab;
    exStruct.waves.opticNerveMask = opticNerveMask;
    exStruct.waves.retinaBoundMask = retinaBoundMask;
else
    exStruct.waves.relDistanceTab = relDistanceTab;
    exStruct.waves.opticNerveMask = exStruct.waves.opticNerveMask;
    exStruct.waves.retinaBoundMask = exStruct.waves.retinaBoundMask;
end

save(exStructPath, "exStruct", '-v7.3');
%% clean up
nImages = ij.WindowManager.getImageCount;
for i = 1: nImages
    stackImpWaves = ij.IJ.getImage();
    stackImpWaves.changes = false;
    stackImpWaves.close;
end
end