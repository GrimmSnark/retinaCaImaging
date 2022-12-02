function prepRetinaCalcium(filePath)
% Function creates dF_F and metadata .mat file for a particular recording

% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps3_100ms-freq.nd2';
% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps1.nd2';
% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_protease_ret5_time-laps4_plus-auto-fluo.nd2';

%% defaults
calciumChan = 1;

%% read in nd2 file
[imStack, metaData] = readFLAMEData(filePath);

%% read in images

% split into channels
imStack = reshape(imStack, size(imStack,1), size(imStack,1), length(metaData.colours.emWavelength), []);

imStackCal = squeeze(imStack(:,:,calciumChan,:));

% downsample
if metaData.image.pixelNum > 512
    for xx = 1:size(imStackCal,3)

        downSampleFactor = metaData.image.pixelNum/512;
        downSampleFactor = 1/downSampleFactor;

        % downsample to 512
        imStackResize(:,:,xx) = imresize(imStackCal(:,:,xx), downSampleFactor);
    end
else
    imStackResize = imStackCal;
end

imStackLinear = reshape(imStackResize, [], size(imStackResize,3));

%% process

%% get the baseline traces
parfor_progress(length(imStackLinear));
parfor i = 1:length(imStackLinear)

    %% fit exp curve to remove bleaching
    y = double(imStackLinear(i,:));
    x = 1:length(y);
    [~, yBaselined(i,:)] = fitexpCurve(x, y);

    %% baseline subtraction
    highpassFilteredTrace(i,:) = baselinePercentileFilter(yBaselined(i,:)', 1/framePeriod ,30);
    parfor_progress;
end
parfor_progress(0);

% create dF/F traces per pixel
dF = (yBaselined-highpassFilteredTrace)./highpassFilteredTrace;
dFRescale = uint16(rescale(dF)* 65536);
dFStack = reshape(dFRescale, 512, 512, []);

%% save stack
[folderParts, name ] = fileparts(filePath);

saveastiff(dFStack, fullfile(folderParts, [name '_dF_F.tif']));
save(fullfile(folderParts, [name '_ExStruct.mat']), 'metaData');
% bfsave(dFStack, fullfile(folderParts, [name(1:end-4) '_dF_F.tif']), 'dimensionOrder', 'XYTCZ');

end