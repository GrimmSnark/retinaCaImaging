function [tifStack,xyShifts, options_nonrigid] = imageRegistration(tifStack,imageRegistrationMethod,spatialResolution,filterCutoff,templateImage, numChunks)
% [tifStack,xyShifts] = imageRegistration(tifStack,imageRegistrationMethod,spatialResolution,filterCutoff,templateImage)
% Registers imaging stack to a template image using either a DFT-based subpixel method ('subMicronMethod') or a
% non rigid NoRMCorr ('nonRigid'). The template image can be directly specified, or the
% the first 100 images are used. If the filter cutoff is empty, then no spatial filtering is done to the images.
%
% by David Whitney (david.whitney@mpfi.org), Max Planck Florida Institute, 2017.
% Modifed by Michael Savage 2022 for proper GPU utilisation


if(nargin<2) || isempty(imageRegistrationMethod), imageRegistrationMethod = 'subMicronMethod'; end % can be either subMicronMethod or downsampleReg
if(nargin<3) || isempty(spatialResolution), spatialResolution = 1.3650; end % in microns per pixel
if(nargin<4) || isempty(filterCutoff), filterCutOff  = [20,150];   end % [lowpass cutoff, highpass cutoff] in units of microns
if(nargin<5), templateImage = [];         end % templateImage is ignored when empty
if(nargin<6) || isempty(filterCutoff), numChunks = 2; end
imgsForTemplate     = [1:100];                % how many images to use for the template
useSpatialFiltering = ~isempty(filterCutOff); % spatially filters the images in an attempt to reduce noise that may impair registration
t0=tic;

options_nonrigid = [];

% Generate a spatially filtered template
if(isempty(templateImage))
    templateImg = uint16(mean(tifStack(:,:,imgsForTemplate),3));
else
    templateImg = templateImage;
end

if(useSpatialFiltering)
    templateImg = real(bandpassFermiFilter(templateImg,-1,filterCutOff(2),spatialResolution));        % Lowpass filtering step
    templateImg = imfilter(templateImg,fspecial('average',round(filterCutOff(1)/spatialResolution))); % Highpass filtering step
end

%% non-rigid NoRMCorr (Flatiron institute)

if strcmp(imageRegistrationMethod, 'nonRigid')
    options_nonrigid = NoRMCorreSetParms('d1',size(tifStack,1),'d2',size(tifStack,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',29,'max_dev',3,'us_fac',50,'init_batch',200);
    [tifStack,xyShifts,~,options_nonrigid,~] = normcorre_batch(tifStack,options_nonrigid, templateImage);
end

%% sub micron DTF rigid correction
if strcmp(imageRegistrationMethod, 'subMicronMethod')
    % Register each image to the template
    numberOfImages = size(tifStack,3);
    xyShifts       = zeros(2,numberOfImages);

    %% GPU version
    if gpuDeviceCount == 1 % tries for GPU  version, which will only work if nvidia CUDA installed
        %     if gpuDeviceCount == 2 % tries for GPU  version, which will only work if nvidia CUDA installed
        % chunk up tifStack if larger than GPU allowance (intmax('int32'))
        t = tic;
        if numel(tifStack) > intmax('int32')

            nblocks = ceil( numel(tifStack)/double(intmax('int32')))+1;
            xyShifts       = [];

            %%

            % get the right chunk numbers
            chunkNo = round(size(tifStack,3)/nblocks);

            % error corrects
            chunkStart = 1:chunkNo:size(tifStack,3);
            chunkEnd = chunkStart + chunkNo-1;
            chunkEnd(chunkEnd >= size(tifStack,3)) = [];
            chunkEnd(end+1) = size(tifStack,3);

            chunkStart(chunkStart >= size(tifStack,3)) = [];
            templateImgGPU = gpuArray(templateImg);

            for cc= 1:nblocks

                % create GPU arrays
                xyShiftsGPU       = gpuArray(zeros(2,chunkEnd(cc)- chunkStart(cc)+1)); % create appropriate size for this chunk
                tifStackGPU = gpuArray(tifStack(:,:,chunkStart(cc): chunkEnd(cc)));

                disp(['Starting to calculate frame shifts using GPU  chunk ' num2str(cc) ' of ' num2str(nblocks)]);

                for ii = 1:size(tifStackGPU,3)
                    % Get current image to register to the template image and pre-process the current frame.

                    sourceImgGPU = tifStackGPU(:,:,ii);
                    if(useSpatialFiltering)
                        sourceImgGPU = real(bandpassFermiFilterGPU(sourceImgGPU,-1,filterCutOff(2),spatialResolution));        % Lowpass filtering step
                        sourceImgGPU = imfilter(sourceImgGPU,fspecial('average',round(filterCutOff(1)/spatialResolution))); % Highpass filtering step
                    end

                    % Determine offsets to shift image
                    [~,output2]=subMicronInPlaneAlignmentGPU(templateImgGPU,sourceImgGPU);
                    xyShiftsGPU(:,ii) = output2(3:4);

                end
                xyShifts = [xyShifts gather(xyShiftsGPU)];
                tifStack(:,:,chunkStart(cc): chunkEnd(cc)) = shiftImageStack(gather(tifStackGPU),gather(xyShiftsGPU([2 1],:)')); % Apply actual shifts to tif stack
            end

        else
            %% full image on GPU
            xyShifts       = zeros(2,numberOfImages);
            templateImgGPU = gpuArray(templateImg);
            tifStackGPU = gpuArray(tifStack);

            disp('Starting to calculate frame shifts using GPU');

%             parfor_progress(numberOfImages);
            parfor ii = 1:numberOfImages
                % Get current image to register to the template image and pre-process the current frame.

%                 parfor_progress; % get progress in parfor loop

                sourceImgGPU = tifStackGPU(:,:,ii);
                if(useSpatialFiltering)
                    sourceImgGPU = real(bandpassFermiFilterGPU(sourceImgGPU,-1,filterCutOff(2),spatialResolution));        % Lowpass filtering step
                    sourceImgGPU = imfilter(sourceImgGPU,fspecial('average',round(filterCutOff(1)/spatialResolution))); % Highpass filtering step
                end

                % Determine offsets to shift image
                [~,output2]=subMicronInPlaneAlignmentGPU(templateImgGPU,sourceImgGPU);
                xyShiftsGPU(:,ii) = output2(3:4);
            end

            tifStack = shiftImageStack(gather(tifStackGPU),gather(xyShiftsGPU([2 1],:)'));
            xyShifts = gather(xyShiftsGPU);

        end
        %% CPU version
    else % tries for CPU version if GPU CUDA not available
        t = tic;
        disp('Starting to calculate frame shifts using CPU');
        xyShifts       = zeros(2,numberOfImages);

%         parfor_progress(numberOfImages)
        parfor ii = 1:numberOfImages
            % Get current image to register to the template image and pre-process the current frame.

%             parfor_progress; % get progress in parfor loop

            sourceImg = tifStack(:,:,ii);
            if(useSpatialFiltering)
                sourceImg = real(bandpassFermiFilter(sourceImg,-1,filterCutOff(2),spatialResolution));        % Lowpass filtering step
                sourceImg = imfilter(sourceImg,fspecial('average',round(filterCutOff(1)/spatialResolution))); % Highpass filtering step
            end

            % Determine offsets to shift image
            [~,output2]=subMicronInPlaneAlignment(templateImg,sourceImg);
            xyShifts(:,ii) = output2(3:4);
        end

        tifStack = shiftImageStack(tifStack,xyShifts([2 1],:)'); % Apply actual shifts to tif stack
    end
end

toc(t);

timeElapsed = toc(t0);
sprintf('Finished registering imaging data - Time elapsed is %4.2f seconds',timeElapsed);

end