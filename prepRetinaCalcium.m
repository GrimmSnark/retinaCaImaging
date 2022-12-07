function prepRetinaCalcium(filePath, motionCorrFlag, motionCorrectionType, createDFPixelMovieFlag)
% Function creates dF_F and metadata .mat file for a particular recording

% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps3_100ms-freq.nd2';
% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps1.nd2';
% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_protease_ret5_time-laps4_plus-auto-fluo.nd2';

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

noOfImagesForAlignment = 50;

%% read in nd2 file
[imStack, metaData] = readFLAMEData(filePath);

%% read in images

% split into channels
imStack = reshape(imStack, size(imStack,1), size(imStack,1), length(metaData.colours.emWavelength), []);

imStackCal = squeeze(imStack(:,:,calciumChan,:));

%% get save folder locations
[folderParts, name ] = fileparts(filePath);

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
% bfsave(dFStack, fullfile(folderParts, [name(1:end-4) '_dF_F.tif']), 'dimensionOrder', 'XYTCZ');

end