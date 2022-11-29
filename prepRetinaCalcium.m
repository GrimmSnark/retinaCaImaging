function prepRetinaCalcium(filePath)

% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps3_100ms-freq.nd2';
% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps1.nd2';
% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_protease_ret5_time-laps4_plus-auto-fluo.nd2';
%% read in nd2 file
imageStruct = bfopen2(filePath);

%% get the meta data
omeMeta = imageStruct{1, 4};

% get channel num
chanNum = omeMeta.getChannelCount(0);

imNo = omeMeta.getPlaneCount(0);

% get frame times
for i = 1:imNo
    frameTimeObj = omeMeta.getPlaneDeltaT(0,i-1);
    frameTime(i) = double(frameTimeObj.value);
end

framePeriod = mean(diff(unique(frameTime)));

%% read in images
imStack = imageStruct{1,1}(:,1);

% split into channels
imStack = cat(3,imStack{:});
imStack = reshape(imStack, size(imStack,1), size(imStack,1), chanNum, []);

imStackCal = squeeze(imStack(:,:,1,:));

% downsample
for xx = 1:size(imStackCal,3)
    imStackResize(:,:,xx) = imresize(imStackCal(:,:,xx), 0.25);
end

imStackLinear = reshape(imStackResize, [], size(imStackResize,3));

frameLen = size(imStackCal,3);

%% process
parfor_progress(length(imStackLinear));
parfor i = 1:length(imStackLinear)

    %% fit exp curve to remove bleaching
    y = double(imStackLinear(i,:));
    x = 1:length(y);

    [~, yBaselined(i,:)] = fitexpCurve(x, y);

    %% baseline subtraction
    %     [~,percentileFiltCutOff] = estimate_percentile_level(double(imStackLinear(i,:))',frameLen,frameLen);
    %     highpassFilteredTrace(i,:) = baselinePercentileFilter(imStackLinear(i,:)', 1/framePeriod ,30,percentileFiltCutOff);
    highpassFilteredTrace(i,:) = baselinePercentileFilter(yBaselined(i,:)', 1/framePeriod ,30);
    parfor_progress;
end
parfor_progress(0);

% create dF/F traces per pixel

% dF = (double(imStackLinear)-highpassFilteredTrace)./highpassFilteredTrace;
dF = (yBaselined-highpassFilteredTrace)./highpassFilteredTrace;

dFRescale = uint16(rescale(dF)* 65536);

dFStack = reshape(dFRescale, 512, 512, []);

%% save stack

[folderParts, name ] = fileparts(filePath);

saveastiff(dFStack, fullfile(folderParts, [name(1:end-4) '_dF_F.tif']));
end