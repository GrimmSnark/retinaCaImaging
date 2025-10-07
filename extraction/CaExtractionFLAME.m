function exStruct = CaExtractionFLAME(exStruct, baselineSubType, channel2Use)
% Function enacts main Ca trace analysis for movie files, applys motion
% correction shifts which were previously calculated, exacts rawF for cells
% and neuropil, does bleaching and baseline correction, calculates dF/F
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs - exStruct (structure containing all experimental data)
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
%          channel2Use: can specify channel to analyse if there are more
%                       than one recorded channel
%                      (OPTIONAL) default = 1 (green channel)
%
% Outputs - exStruct (updated structure)

%% set defaults

if nargin <2 || isempty(baselineSubType)
    baselineSubType = 1;
end

if nargin <3 || isempty(channel2Use)
    channel2Use = 1; % sets default channel to use if in multi channel recording
end

%% Basic setup of tif stack

% sets up ROI manager for this function
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();

%% read in nd2 file

[~, ~, ext ] = fileparts(exStruct.filePath);

%% read in image file
if strcmp(ext, '.nd2')
    imStack = readFLAMEData(exStruct.filePath);
    %% read in images
    % split into channels
    imStack = reshape(imStack, size(imStack,1), size(imStack,1), length(exStruct.colours.emWavelength), []);
    vol = squeeze(imStack(:,:,channel2Use,:));
    imStack = [];

else
    try
        vol = read_Tiffs(exStruct.filePath);
    catch
        vol = readMultipageTifFiles(exStruct.filePath);
    end
end


% apply imageregistration shifts if there are shifts to apply
if isprop(exStruct, 'options_nonrigid') && ~isempty(exStruct.options_nonrigid) % if using non rigid correctionn
    registeredVol = apply_shifts(vol,exStruct.xyShifts,exStruct.options_nonrigid);
    % elseif  ~isempty(exStruct.xyShifts) && isfield(exStruct, 'xyShifts')
elseif isfield(exStruct, 'xyShifts')  && ~isempty(exStruct.xyShifts)
    registeredVol = shiftImageStack(vol,exStruct.xyShifts([2 1],:)'); % Apply actual shifts to tif stack
else % if there are no motion correction options, ie the image stack loaded is already motion corrected
    registeredVol = vol;
end

% transfers to FIJI
registeredVolMIJI = MIJ.createImage( 'Registered Volume', registeredVol,true);

%% start running raw trace extraction

% allocate fields
exStruct.cells.rawF = [];
exStruct.cells.rawF_neuropil =[];
exStruct.cells.xPos = zeros(exStruct.cells.cellCount,1);
exStruct.cells.yPos = zeros(exStruct.cells.cellCount,1);

for x = 1:exStruct.cells.cellCount
    % Select cell ROI in ImageJ/FIJI
    fprintf('Processing Cell %d\n',x)

    % Get cell ROI name and parse out (X,Y) coordinates
    RC.select(x-1); % Select current cell
    currentROI = RC.getRoi(x-1);

    exStruct.cells.yPos(x) = currentROI.getYBase;
    exStruct.cells.xPos(x) = currentROI.getXBase;


    % Get the fluorescence timecourse for the cell and neuropil ROI by
    % using ImageJ's "z-axis profile" function.
    for isNeuropilROI = 0:1
        ij.IJ.getInstance().toFront();

        plotTrace = ij.plugin.ZAxisProfiler.getPlot(registeredVolMIJI);
        RT(:,1) = plotTrace.getXValues();
        RT(:,2) = plotTrace.getYValues();

        if isNeuropilROI
            %RC.setName(sprintf('Neuropil ROI %d',i));
            exStruct.cells.rawF_neuropil(x,:) = RT(:,2);
        else
            %RC.setName(sprintf('Cell ROI %d',i));
            exStruct.cells.rawF(x,:) = RT(:,2);
            RC.select((x-1)+exStruct.cells.cellCount); % Now select the associated neuropil ROI
        end
    end
end

%% subtract neuropil signal
exStruct.cells.baseline  = 0*exStruct.cells.rawF;

switch baselineSubType

    case 1
        %% expotential bleaching fit then rolling ball baseline median percentile filter

        % Define a moving and fluorescence baseline using a percentile filter
        for q =1:exStruct.cells.cellCount
            fprintf('Computing baseline for cell %d using Exp bleach fit and baseline subtraction\n', q);

            % fit exp curve to remove bleaching
            [~, exStruct.cells.baslineBleaching(q,:)] = fitexpCurve(1:size(exStruct.cells.rawF,2), exStruct.cells.rawF(q,:));

            % baseline subtraction
            try
                exStruct.cells.baseline(q,:) = baselinePercentileFilter(exStruct.cells.baslineBleaching(q,:)', exStruct.cells.fps ,30);
            catch
                exStruct.cells.baseline(q,:) = baselinePercentileFilter(exStruct.cells.baslineBleaching(q,:)', exStruct.rate ,30);
            end

        end

        % computer delta F/F traces
        exStruct.cells.dF = (exStruct.cells.baslineBleaching-exStruct.cells.baseline)./exStruct.cells.baseline;

    case 2
        %% Annulus neuropil subtraction and kernal density estimation for percentile filter

        % Compute the neuropil-contributed signal from our cells, and eliminate
        % them from our raw cellular trace

        [exStruct.cells.correctedF, exStruct.cells.neuropCorrPars]=estimateNeuropil(exStruct.cells.rawF,exStruct.cells.rawF_neuropil); % neuropil subtraction calculated from Dipoppa et al 2018

        % Define a moving and fluorescence baseline using a percentile filter

        for q =1:exStruct.cells.cellCount
            fprintf('Computing baseline for cell %d \n', q);

            % Compute a moving baseline with a 60s percentile lowpass filter smoothed by a 60s Butterworth filter
            offset = 5000; % helps get around dividing my zero errors...
            [~,percentileFiltCutOff(q)] = estimate_percentile_level((exStruct.cells.correctedF(q,:)'+offset),size(registeredVol,3),size(registeredVol,3));

            lowPassFiltCutOff    = 30; %in seconds

            try
                exStruct.cells.baseline(q,:)  = baselinePercentileFilter((exStruct.cells.correctedF(q,:)'+offset),exStruct.rate,lowPassFiltCutOff,percentileFiltCutOff(q));
            catch
                 exStruct.cells.baseline(q,:)  = baselinePercentileFilter((exStruct.cells.correctedF(q,:)'+offset),exStruct.fps,lowPassFiltCutOff,percentileFiltCutOff(q));
            end

        end

        % store percentile rank filter cutoff for each cell
        exStruct.cells.percentileFiltCutOff = percentileFiltCutOff;

        % store corrected baseline (i.e remove offset)
        exStruct.cells.baselineCorrected = exStruct.cells.baseline-offset;

        % computer delta F/F traces
        exStruct.cells.dF = ((exStruct.cells.correctedF+offset)-exStruct.cells.baseline)./exStruct.cells.baseline;
end
end