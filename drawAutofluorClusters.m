function drawAutofluorClusters(clusterImg, exStructPath)

%% open FIJI
% initalize MIJI 
intializeMIJ;
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();

%% load in exStruct
exStruct = load(exStructPath);

try
    exStruct = exStruct.exStruct;
catch
    exStruct = exStruct.metaData;
end

%% load image

MIJ.run('Open...', ['path=' clusterImg]);
stackImp = ij.IJ.getImage();
MIJ.run("Split Channels");

imageNamesJava = MIJ.getListImages;
textChanGrand =[];

for x = 1:length(imageNamesJava)
windowNames{x} = char(imageNamesJava(x));
MIJ.selectWindow(windowNames{x});
MIJ.run('Gaussian Blur...', 'sigma=3');
MIJ.run("Subtract Background...", "rolling=10 sliding");
MIJ.run("Enhance Contrast...", "saturated=0.35");

textChan = sprintf('c%i=%s ',x, windowNames{x});
textChanGrand = [textChanGrand textChan];
end

MIJ.run("Merge Channels...", [textChanGrand 'create']);
ij.IJ.runMacro('setTool("polygon");');

%% choose ROI shapes

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
clusterCellMask = createLabeledROIFromImageJPixels([pixelIm pixelIm] ,ROIobjects);

shapeProps = regionprops("table", clusterCellMask,"Area", "BoundingBox","Centroid", "SubarrayIdx","PixelIdxList", "ConvexHull");

%% save into exStruct
exStruct.clusterCells.mask = clusterCellMask;
exStruct.clusterCells.clusterProps = shapeProps;

%% Clean up windows
stackImp = ij.IJ.getImage();
stackImp.changes = false;
stackImp.close;
MIJ.closeAllWindows;

end