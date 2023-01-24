function CaAnalysisFLAME(exStructPath, baselineSubType)
% function which runs main analysis on calcium imaging data recorded on
% FLAME FN1 system. This function requires user input defined cell
% ROIs and calculates dF/F with baseline subtraction.
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs-  exStructPath: filepath for exStruct.mat to process
%
%          baselineSubType - Switch case for calcium imaging baseline
%                            subtraction type DEFAULT == 1
%                            1: Expotential bleaching fit then rolling ball
%                               baseline median percentile filter (use for
%                               highly packed cells, ie retina)
%                            2: Annulus neuropil subtraction and kernal 
%                               density estimation for percentile filter (
%                               use for loosely packed cells, ie culture or
%                               in-vivo brain)
%
% Output- exStruct: saves structure containing all experiment info

%% set defaults
if nargin <2 || isempty(baselineSubType)
    baselineSubType = 1;
end

%% open FIJI
% initalize MIJI and get ROI manager open
intializeMIJ;
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();
RC.reset();

%% load metaStruct
exStruct = load(exStructPath);
exStruct = exStruct.exStruct;

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

exStruct.cells.labeledCellROI = createLabeledROIFromImageJPixels([exStruct.image.pixelNum exStruct.image.pixelNum], cellROIs);
exStruct.cells.labeledNeuropilROI = createLabeledROIFromImageJPixels([exStruct.image.pixelNum exStruct.image.pixelNum], neuropilROIs);
exStruct.cells.averageROIRadius = averageROIRadius;

% does calcium trace extraction
exStruct = CaExtractionFLAME(exStruct, baselineSubType);


% save data
save(exStructPath, "exStruct", '-v7.3');

% Clean up windows
MIJ.closeAllWindows;

end