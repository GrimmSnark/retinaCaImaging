function CaAnalysisFLAME(exStructPath)
% function which runs main analysis on calcium imaging data recorded on
% FLAME FN1 system. This function requires user input defined cell
% ROIs and calculates dF/F with baseline subtraction.
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs-  exStructPath: filepath for exStruct.mat to process
%
% Output- exStruct: saves structure containing all experiment info

%% open FIJI
% initalize MIJI and get ROI manager open
intializeMIJ;
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();

%% load metaStruct
exStruct = load(exStructPath);
exStruct = exStruct.metaData;

%% get filepath root
filePathRoot = exStructPath(1:end-13);

%% start analysis

% load in pointers to ROI manager
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();


% load in ROI file
if exist([filePathRoot '.zip'])
    RC.runCommand('Open', [filePathRoot '.zip']); % opens zip file
else
    disp(['No ROI file found in "' filePathRoot '.zip" Please run selectROIs.m']);
    return
end
ROInumber = RC.getCount();
disp(['You have selected ' num2str(ROInumber) ' ROIs, moving onto analysis']);

% calculate average cell ROI radius
averageROIRadius = calculateNeuropilRoiRadius(RC.getRoisAsArray);

% get cell ROI radius, match neuropil ROI radius
generateNeuropilROIs(RC.getRoisAsArray,(averageROIRadius*4)); % generates neuropil surround ROIs


%% load in experimentStructure and begin trace extraction

% feeds in data into experiement structure
exStruct.cellCount = ROInumber;
ROIobjects = RC.getRoisAsArray;
cellROIs = ROIobjects(1:ROInumber);
neuropilROIs = ROIobjects(ROInumber+1:end);

exStruct.labeledCellROI = createLabeledROIFromImageJPixels([exStruct.image.pixelNum exStruct.image.pixelNum], cellROIs);
exStruct.labeledNeuropilROI = createLabeledROIFromImageJPixels([exStruct.image.pixelNum exStruct.image.pixelNum], neuropilROIs);
exStruct.averageROIRadius = averageROIRadius;

% does calcium trace extraction
exStruct = CaExtractionFLAME(exStruct);


% save data
save(exStructPath, "exStruct", '-v7.3');

% Clean up windows
MIJ.closeAllWindows;

end