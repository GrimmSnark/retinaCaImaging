function prepRetinaCalcium(filePath, motionCorrFlag, motionCorrectionType, createDFPixelMovieFlag)
% This function does basic preprocessing of t series imaging data including
% meta data, image registration,  creates dF_F pixelwise movie and summary
% SD images from .nd2 files from FLAME system
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Input- filePath: image data .nd2 file from FLAME system
%
%        motionCorrFlag- Flag do motion correction on raw images
%
%        motionCorrectionType - DFT-based subpixel method
%                               ('subMicronMethod')
%                             - non-rigid NoRM Corr registration
%                               ('nonRigid')
%
%        createDFPixelMovieFlag: flag for creating pixel wise DF_F stack
%        and takes up time/space 0 = not saved, 1 = saved (DEFAULT)

%% defaults
calciumChan = 1;

if nargin < 2 || isempty(motionCorrFlag)
    motionCorrFlag = 1;
end

if nargin < 3 || isempty(motionCorrectionType)
    motionCorrectionType = 'subMicronMethod';
end


if nargin < 4 || isempty(createDFPixelMovieFlag)
    createDFPixelMovieFlag = 1;
end

noOfImagesForAlignment = 50; % number of brightest images used for motion correction template image

%% get save folder locations
[folderParts, name, ext ] = fileparts(filePath);

%% read in image file
if strcmp(ext, '.nd2')
    [imStack, metaData] = readFLAMEData(filePath);
    % split into channels
    imStack = reshape(imStack, size(imStack,1), size(imStack,1), length(metaData.colours.emWavelength), []);
    imStackCal = squeeze(imStack(:,:,calciumChan,:));
else
    imStackCal = read_Tiffs(filePath); 
    metaData = [];
    metaData.filePath = filePath;
    metaData.rate = 10;
    metaData.image.pixelNum = size(imStackCal, 1);
end

%% motion correction

if motionCorrFlag == 1
    % we want to use the brightness average
    % get image brightness in stack
    imageBrightness = squeeze(mean(imStackCal,[1,2]));

    % sort image brightness
    [~, imageBrightnessIndx] = sort(imageBrightness);
    % make average image based on noOfImages brightest images
    templateImageForReg = uint16(mean(imStackCal(:,:,imageBrightnessIndx(1:noOfImagesForAlignment-1)),3));

    disp(['Starting image registration using ' motionCorrectionType ]);
    [imStackCal,xyShifts, options_nonrigid] = imageRegistration(imStackCal,motionCorrectionType, metaData.image.pixelSize, [], templateImageForReg);
    metaData.xyShifts = xyShifts;
    metaData.options_nonrigid = options_nonrigid;
else
    metaData.xyShifts = [];
    metaData.options_nonrigid = [];
end

%% create DF/F image stack

if createDFPixelMovieFlag == 1
    [dFStack, metaData] = createDFPixelImageStack(imStackCal, metaData);

    % create SD of DF movie
    imageDF_SD = std(double(dFStack),[],3);
    imageDF_SD = uint16(mat2gray(imageDF_SD) * 65535);
    saveastiff(imageDF_SD, fullfile(folderParts, [name '_dF_SD.tif']));
end
%% save meta and SD image

imageSD = std(double(imStackCal),[],3);
imageSD = uint16(mat2gray(imageSD) * 65535);
saveastiff(imageSD, fullfile(folderParts, [name '_SD.tif']));

save(fullfile(folderParts, [name '_ExStruct.mat']), 'metaData');
end