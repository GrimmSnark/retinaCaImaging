function prepRetinaCalcium()

% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps3_100ms-freq.nd2';
% filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps1.nd2';
% %% read in nd2 file
% imageStruct = bfopen(filePath);
% imStack = imageStruct{1,1}(:,1);
% imStackCal = cat(3,imStack{:});


filePath = 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_protease_ret5_time-laps4_plus-auto-fluo.nd2';
%% read in nd2 file
imageStruct = bfopen(filePath);
imStack = imageStruct{1,1}(:,1);
imStack = cat(3,imStack{:});

imStackCal = imStack(:,:,1:2:end);

% downsample
for xx = 1:size(imStackCal,3)
imStackResize(:,:,xx) = imresize(imStackCal(:,:,xx), 0.25);
end

% imStackLinear = reshape(imStack, [], size(imStack,3));
imStackLinear = reshape(imStackResize, [], size(imStackResize,3));

frameLen = size(imStackCal,3);

 parfor_progress(length(imStackLinear));
% timeStart = tic;
% for i = 1:length(imStackLinear)
% 
%     %% baseline subtraction
%     highpassFilteredTrace(i,:) = baselinePercentileFilter(imStackLinear(i,:)',3,30);
%     prcdone(i,length(imStackLinear),'loop time' ,1,timeStart);
% %     parfor_progress;
% end

parfor i = 1:length(imStackLinear)

    %% baseline subtraction
    [~,percentileFiltCutOff] = estimate_percentile_level(double(imStackLinear(i,:))',frameLen,frameLen);
    highpassFilteredTrace(i,:) = baselinePercentileFilter(imStackLinear(i,:)',3,30,percentileFiltCutOff);
    parfor_progress;
end
parfor_progress(0);

dF = (double(imStackLinear)-highpassFilteredTrace)./highpassFilteredTrace;
dFRescale = uint16(rescale(dF)* 65536);

dFStack = reshape(dFRescale, 512, 512, []);

saveastiff(dFStack, 'C:\Data\mouse\calcium\FN1\09-11_Cal_520_ret3_time-laps1_dF_F.tif' )
%% SD projection for ROI

end