function drawAutofluorClusters(clusterImg, exStructPath, imgProcess)


if nargin < 3 || isempty(imgProcess)
    imgProcess = 1;
end

%% open FIJI
% initalize MIJI
intializeMIJ;

RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();
RC.reset();

%% load in exStruct

if nargin <2 || isempty(exStructPath)
    exStructPath = [clusterImg(1:end-4) '_ExStruct.mat'];
end

exStruct = load(exStructPath);

try
    exStruct = exStruct.exStruct;
catch
    exStruct = exStruct.metaData;
end

%% load image

MIJ.run('Open...', ['path=' clusterImg]);
stackImp = ij.IJ.getImage();

stackDim = stackImp.getDimensions();

% if it is a multi-colour image
if stackDim(3) > 1
    MIJ.run("Split Channels");

    imageNamesJava = MIJ.getListImages;
    textChanGrand =[];

    for x = 1:length(imageNamesJava)

        windowNames{x} = char(imageNamesJava(x));

        if stackImp.getNFrames > 1
            ij.IJ.selectWindow([windowNames{x}]);

            MIJ.run("Z Project...", "projection=[Max Intensity]");

            if imgProcess == 1
                MIJ.run('Gaussian Blur...', 'sigma=3');
                MIJ.run("Subtract Background...", "rolling=10 sliding");
                MIJ.run("Enhance Contrast...", "saturated=0.35");
            end

            textChan = sprintf('c%i=MAX_%s ',x, windowNames{x});
        else
            MIJ.selectWindow([windowNames{x}]);

            if imgProcess == 1
                MIJ.run('Gaussian Blur...', 'sigma=3');
                MIJ.run("Subtract Background...", "rolling=10 sliding");
                MIJ.run("Enhance Contrast...", "saturated=0.35");
            end

            textChan = sprintf('c%i=%s ',x, windowNames{x});
        end

        textChanGrand = [textChanGrand textChan];
    end

    MIJ.run("Merge Channels...", [textChanGrand 'create']);

end
ij.IJ.runMacro('setTool("polygon");');

%% choose ROI shapes for cluster cells

% Sets up diolg box to allow for user input to choose cell ROIs
opts.Default = 'Continue';
opts.Interpreter = 'tex';

questText = [{'Draw ROIs around the cluster cell locations'} ...
    {'Add ROI to ROI Manager by Ctrl + T'} {''} ...
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

pixelIm = stackImp.getWidth;

% full size image
clusterCellMask = createLabeledROIFromImageJPixels([pixelIm pixelIm] ,ROIobjects);
shapeProps = regionprops("table", clusterCellMask,"Area", "BoundingBox","Centroid", "SubarrayIdx","PixelIdxList", "ConvexHull");

% downsampled image
rz = 512/pixelIm;
clusterCellMaskDownsize = imresize(clusterCellMask,rz, 'nearest');
shapePropsDownsize = regionprops("table", clusterCellMaskDownsize,"Area", "BoundingBox","Centroid", "SubarrayIdx","PixelIdxList", "ConvexHull");

%% choose ROI shapes for blood vessel plexus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RC.reset();
% Sets up diolg box to allow for user input to choose cell ROIs
opts.Default = 'Continue';
opts.Interpreter = 'tex';

questText = [{'Draw ROI around the blood vessel plexus'} ...
    {'Add ROI to ROI Manager by Ctrl + T'} {''} ...
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

pixelIm = stackImp.getWidth;

% full size image
bloodVesselMask = createLabeledROIFromImageJPixels([pixelIm pixelIm] ,ROIobjects);
shapePropsBloodVessel = regionprops("table", bloodVesselMask,"Area", "BoundingBox","Centroid", "SubarrayIdx","PixelIdxList", "ConvexHull");

% downsampled image
rz = 512/pixelIm;
bloodVesselMaskDownsize = imresize(bloodVesselMask,rz, 'nearest');
shapePropsBloodVesselDownsize = regionprops("table", bloodVesselMaskDownsize,"Area", "BoundingBox","Centroid", "SubarrayIdx","PixelIdxList", "ConvexHull");





%% save into exStruct
exStruct.clusterCells.mask = clusterCellMask;
exStruct.clusterCells.clusterProps = shapeProps;

exStruct.clusterCells.maskDS = clusterCellMaskDownsize;
exStruct.clusterCells.clusterPropsDS = shapePropsDownsize;

exStruct.bloodVesselPlexus.mask = bloodVesselMask;
exStruct.bloodVesselPlexus.bloodVesselProps = shapePropsBloodVessel;

exStruct.bloodVesselPlexus.maskDS = bloodVesselMaskDownsize;
exStruct.bloodVesselPlexus.bloodVesselPropsDS = shapePropsBloodVesselDownsize;


% save data
save(exStructPath, "exStruct", '-v7.3');

%% Clean up windows

nImages = ij.WindowManager.getImageCount;

for n = 1:nImages
    stackImp = ij.IJ.getImage();
    stackImp.changes = false;
    stackImp.close;
end

% ij.WindowManager.closeAllWindows;

end