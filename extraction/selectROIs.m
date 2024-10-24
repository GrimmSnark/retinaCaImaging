function selectROIs(exStructPath)
% Function loads in any standard deviation images for a recording and plugs
% them into FIJI to allow user defined cell ROIs. Saves .zip file ROIs for
% Ca extraction
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs-  exStructPath: filepath for exStruct.mat to process

%% set default

killFlag = 0;

if nargin < 1 || isempty(exStructPath)
   [file, path] = uigetfile({'*.mat'},...
                          'Image File Selector');

   exStructPath = fullfile(path,file);
end

% get the appropriate magnification for ROI image
screenDim = get(0,'ScreenSize');
if screenDim(3) > 2000
    magSize = 300; % magnification for image viewing in %
else
    magSize = 200; % magnification for image viewing in %
end

%% open FIJI
% initalize MIJI and get ROI manager open
intializeMIJ;
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();

%% load metaStruct
exStruct = load(exStructPath);

try
    exStruct = exStruct.metaData;
catch
    exStruct = exStruct.exStruct;
end

%% get filepath root
filePathRoot = exStructPath(1:end-13);

%% Create the appropriate images for ROI extraction

% finds all the relevant images for ROI choosing
files = dir([filePathRoot '*SD*.tif']);

if size(files,1) ==1 % if single channel recording
    imageROI = read_Tiffs(fullfile(files(1).folder, files(1).name),1); % reads in average image
    imageROI = imadjust(imageROI); % saturate image to make neural net prediction better
end

% get both SD images if you ran the dF/F pixelwise movie creation
if size(files,1) >1 % if multiple channel recording

    % read all SD images in for the recording
    for xx = 1:length(files)
        eval(['imageROI' num2str(xx) ' = read_Tiffs(''' fullfile(files(xx).folder, files(xx).name) ''',1);']);
    end

    % get the image sizes
     for xx = 1:length(files)
        eval(['imSize(' num2str(xx) ') = size(imageROI' num2str(xx) ',1);']);
    end
    
    % get smaller and find the magnfication factor
    [~, minIn] = min(imSize);
    [~, maxIn] = max(imSize);
    magFactor = imSize(maxIn)/imSize(minIn);

    % resize image
    eval(['imageROI' num2str(minIn) ' = imresize(imageROI' num2str(minIn) ',' num2str(magFactor) ');']);

    % concatenate images
    imageROI = cat(3, imadjust(imageROI1), imadjust(imageROI2));
end

%% open image in FIJI
MIJ.createImage(imageROI);

pause(0.2);
ij.IJ.run('Set... ', ['zoom=' num2str(magSize) ' x=10 y=50']);
ij.IJ.run('Enhance Contrast', 'saturated=0.35');

% open ROI tool
MIJ.run("Cell Magic Wand Tool");

%% Deal with ROI selection

% Check if there are already ROIs selected for this recording
ROI_File = dir([filePathRoot '*.zip']);

if ~isempty(ROI_File)
    disp([[fullfile(ROI_File.folder, ROI_File.name)]  ' contains a valid ROI file!']);
    RC.runCommand('Open', [fullfile(ROI_File.folder, ROI_File.name)]); % opens ROI file
    
    % Query user if you want to use previously chosen ROIs
    answer = MFquestdlg([0.5,0.5], 'Would you like to use previously chosen ROIs?', ...
        'Choose your ROIs', ...
        'Yes','No', 'Yes');
    % Handle response
    switch answer
        case 'Yes'
            
        case 'No'
            RC.runCommand('Delete'); % resets ROIs if you select clear all
        case ''

    end
end

% Sets up diolg box to allow for user input to choose cell ROIs
opts.Default = 'Continue';
opts.Interpreter = 'tex';

questText = [{'Check ROIs and remove any covering multiple cells'} ...
    {'Select the Cell Magic Wand and hold shift and click to add to ROI manager'} {''} ...
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
        RC.runCommand('Save', [filePathRoot '.zip']); % saves zip file
        
    case 'Exit Script' % if you want to exit and end
        killFlag = 1;
        MIJ.closeAllWindows;
        return
    case ''
        killFlag = 1; % if you want to exit and end
        MIJ.closeAllWindows;
        return
end

%% clean up
nImages = ij.WindowManager.getImageCount;
for i = 1: nImages
    stackImpWaves = ij.IJ.getImage();
    stackImpWaves.changes = false;
    stackImpWaves.close;
end

RC.close;
%% update and save exStruct

exStruct.cells.cellCount = ROInumber;
save(exStructPath, 'exStruct');
end