function [dFStack, imageMetaData] = createDFPixelImageStack(imgStack, imageMetaData)
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
%          imageMetaData - metaData structure with the original DF/F image stack
%                     values as double array

%% downsample  to make computation easier

intializeMIJ;

if imageMetaData.image.pixelNum > 512
    for xx = 1:size(imgStack,3)

        downSampleFactor = imageMetaData.image.pixelNum/512;
        downSampleFactor = 1/downSampleFactor;

        % downsample to 512
        imStackResize(:,:,xx) = imresize(imgStack(:,:,xx), downSampleFactor);
        imageMetaData.downsampledRes = downSampleFactor * imageMetaData.image.pixelSize;
    end
else
    imStackResize = imgStack;
end

imp = MIJ.createImage('stack', imStackResize, 0);
bc = emblcmci.BleachCorrection;
bc.setCorrectionMethod(1);
bc.setHeadlessProcessing(1);
impcorrected = bc.doCorrection(imp);

for i =1:size(imStackResize,3)
    impcorrected.setSlice(i-1);
    imStackBleachCorr(:,i) =impcorrected.getProcessor().getPixels;
end

imStackBleachCorr = reshape(imStackBleachCorr, size(imStackResize));
imStackBleachCorr = permute(imStackBleachCorr, [2 1 3]);

% imStackLinear = reshape(imStackResize, [], size(imStackResize,3));
imStackLinear = reshape(imStackBleachCorr, [], size(imStackResize,3));

%% process

%% get the baseline traces
framePeriod = imageMetaData.framePeriod;

if gpuDeviceCount == 1

    %% run computation on GPU
    imStackLinearGPU = gpuArray(imStackLinear);
    %     parfor_progress(length(imStackLinear));
    tic
    disp('Working on pixelwise dF/F image, will take some time...')
    parfor i = 1:length(imStackLinear)
        %   for i = 1:length(imStackLinear)


        %% fit exp curve to remove bleaching on GPU
        y = double(imStackLinearGPU(i,:));
        x = 1:length(y);

        % added to avoid some weird error in how the baselining
        % messes up the traces, take the original instead of the
        % corrected
        try
            [~, yBaselined(i,:)] = fitExpCurveGPU(x, y);
        catch
            [r, c] = ind2sub(size(imStackResize, [1 2]), i);
            [~, yBaselined(i,:)] = fitExpCurveGPU(x, gpuArray(double(squeeze(imStackResize(r,c,:))')));
        end

        %% baseline subtraction
        highpassFilteredTrace(i,:) = baselinePercentileFilter2(yBaselined(i,:)', 1/framePeriod ,30);
        %         highpassFilteredTrace(i,:) = baselinePercentileFilter2(double(imStackLinearGPU(i,:))', 1/framePeriod ,30);


        %         prcdone(i,length(imStackLinear),'Baseline on GPU',10 ,tStart);
        %         parfor_progress;
    end
    %     parfor_progress(0);

    highpassFilteredTrace = gather(highpassFilteredTrace);
    yBaselined = gather(yBaselined);

    toc;
else
    % run computation on CPU
    %     parfor_progress(length(imStackLinear));
    parfor i = 1:length(imStackLinear)
        % fit exp curve to remove bleaching
        y = double(imStackLinear(i,:));
        x = 1:length(y);
        [~, yBaselined(i,:)] = fitexpCurve(x, y);

        %% baseline subtraction
        highpassFilteredTrace(i,:) = baselinePercentileFilter(yBaselined(i,:)', 1/framePeriod ,30);
        %         parfor_progress;
    end
    %     parfor_progress(0);
end

disp('Completed pixelwise dF/F image')


% close parpool to save space on RAM...
poolobj = gcp('nocreate');
delete(poolobj);

% create dF/F traces per pixel
% dF = (double(imStackLinear)-highpassFilteredTrace)./highpassFilteredTrace;
dF = (yBaselined-highpassFilteredTrace)./highpassFilteredTrace;
dFRescale = uint16(rescale(dF)* 65536);
dFStack = reshape(dFRescale, 512, 512, []);

imageMetaData.dF = reshape(dF, 512, 512, []); % save the actual dF/F signal for further analysis

end