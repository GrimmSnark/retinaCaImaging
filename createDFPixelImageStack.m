function [dFStack, metaData] = createDFPixelImageStack(imgStack, metaData)
% creates a pixel-wise DF/F time series stack and adds DF/F matrix to
% metaData structure
%
% NB Will take a long time to run as it runs through all pixels in parfor.
% Can be sped up with GPU but not currently implemented
%
% Inputs: imgStack - XYT image stack matrix to use for DF/F calculations
%
%         metaData - metadata structure for the recordings
%
% Outputs: dFStack - uint16 XYT image stack for display purposes (is
%                    rescaled so not compariable for intensity values)
%
%          metaData - metaData structure with the original DF/F image stack
%                     values as double array

%% downsample  to make computation easier
if metaData.image.pixelNum > 512
    for xx = 1:size(imgStack,3)

        downSampleFactor = metaData.image.pixelNum/512;
        downSampleFactor = 1/downSampleFactor;

        % downsample to 512
        imStackResize(:,:,xx) = imresize(imgStack(:,:,xx), downSampleFactor);
        metaData.downsampledRes = downSampleFactor * exStruct.image.pixelSize;
    end
else
    imStackResize = imgStack;
end

imStackLinear = reshape(imStackResize, [], size(imStackResize,3));

%% process

%% get the baseline traces
framePeriod = metaData.framePeriod;

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

metaData.dF = reshape(dF, 512, 512, []); % save the actual dF/F signal for further analysis

end