function exStruct = CaExtractionFLAME(exStruct, channel2Use)
% Function enacts main Ca trace analysis for movie files, applys motion
% correction shifts which were previously calculated, exacts rawF for cells
% and neuropil, does bleaching and baseline correction, calculates dF/F
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs - exStruct (structure containing all experimental data)
%
%          channel2Use: can specify channel to analyse if there are more
%                       than one recorded channel
%                      (OPTIONAL) default = 1 (green channel)
%
% Outputs - exStruct (updated structure)

%% set defaults

if nargin <2 || isempty(channel2Use)
    channel2Use = 1; % sets default channel to use if in multi channel recording
end

%% Basic setup of tif stack

% sets up ROI manager for this function
RM = ij.plugin.frame.RoiManager();
RC = RM.getInstance();

%% read in nd2 file
imStack = readFLAMEData(exStruct.filePath);

%% read in images

% split into channels
imStack = reshape(imStack, size(imStack,1), size(imStack,1), length(exStruct.colours.emWavelength), []);
vol = squeeze(imStack(:,:,channel2Use,:));


% apply imageregistration shifts if there are shifts to apply
if isprop(exStruct, 'options_nonrigid') && ~isempty(exStruct.options_nonrigid) % if using non rigid correctionn
    registeredVol = apply_shifts(vol,exStruct.xyShifts,exStruct.options_nonrigid);
% elseif  ~isempty(exStruct.xyShifts) && isfield(exStruct, 'xyShifts')
elseif isfield(exStruct, 'xyShifts')
    registeredVol = shiftImageStack(vol,exStruct.xyShifts([2 1],:)'); % Apply actual shifts to tif stack
else % if there are no motion correction options, ie the image stack loaded is already motion corrected
    registeredVol = vol;
end

% transfers to FIJI
registeredVolMIJI = MIJ.createImage( 'Registered Volume', registeredVol,true);

%% start running raw trace extraction

% allocate fields
exStruct.rawF = [];
exStruct.rawF_neuropil =[];
exStruct.xPos = zeros(exStruct.cellCount,1);
exStruct.yPos = zeros(exStruct.cellCount,1);

for x = 1:exStruct.cellCount
    % Select cell ROI in ImageJ/FIJI
    fprintf('Processing Cell %d\n',x)
    
    % Get cell ROI name and parse out (X,Y) coordinates
    RC.select(x-1); % Select current cell
    currentROI = RC.getRoi(x-1);
%     [tempLoc1,tempLoc2] = strtok(char(RC.getName(x-1)),'-');
%     experimentStructure.yPos(x) =  str2double(tempLoc1);
%     experimentStructure.xPos(x) = -str2double(tempLoc2);

    exStruct.yPos(x) = currentROI.getYBase;
    exStruct.xPos(x) = currentROI.getXBase;

    
    % Get the fluorescence timecourse for the cell and neuropil ROI by
    % using ImageJ's "z-axis profile" function.
    for isNeuropilROI = 0:1
        ij.IJ.getInstance().toFront();
        
        plotTrace = ij.plugin.ZAxisProfiler.getPlot(registeredVolMIJI);
        RT(:,1) = plotTrace.getXValues();
        RT(:,2) = plotTrace.getYValues();
        
        if isNeuropilROI
            %RC.setName(sprintf('Neuropil ROI %d',i));
            exStruct.rawF_neuropil(x,:) = RT(:,2);
        else
            %RC.setName(sprintf('Cell ROI %d',i));
            exStruct.rawF(x,:) = RT(:,2);
            RC.select((x-1)+exStruct.cellCount); % Now select the associated neuropil ROI
        end
    end
end

%% subtract neuropil signal

% Compute the neuropil-contributed signal from our cells, and eliminate
% them from our raw cellular trace

[exStruct.correctedF_Neuropil, exStruct.neuropCorrPars]=estimateNeuropil(exStruct.rawF,exStruct.rawF_neuropil); % neuropil subtraction calculated from Dipoppa et al 2018

%% Define a moving and fluorescence baseline using a percentile filter
% exStruct.rate      = 1/exStruct.framePeriod; % frames per second
exStruct.baseline  = 0*exStruct.rawF;


for q =1:exStruct.cellCount
    fprintf('Computing baseline for cell %d \n', q);

      %% fit exp curve to remove bleaching
    [~, exStruct.baslineBleaching(q,:)] = fitexpCurve(1:size(exStruct.rawF,2), exStruct.rawF(q,:));

    %% baseline subtraction
     exStruct.baseline(q,:) = baselinePercentileFilter(exStruct.baslineBleaching(q,:)', exStruct.fps ,30);


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Compute a moving baseline with a 60s percentile lowpass filter smoothed by a 60s Butterworth filter
%     offset = 5000; % helps get around dividing my zero errors...
%     [~,percentileFiltCutOff(q)] = estimate_percentile_level((exStruct.correctedF(q,:)'+offset),size(registeredVol,3),size(registeredVol,3));
%     
%    lowPassFiltCutOff    = 30; %in seconds
%     exStruct.baseline(q,:)  = baselinePercentileFilter((exStruct.correctedF(q,:)'+offset),exStruct.rate,lowPassFiltCutOff,percentileFiltCutOff(q));

end

% store percentile rank filter cutoff for each cell
% exStruct.percentileFiltCutOff = percentileFiltCutOff;

% store corrected baseline (i.e remove offset)
% exStruct.baselineCorrected = exStruct.baseline-offset;

%% computer delta F/F traces
% exStruct.dF = ((exStruct.correctedF+offset)-exStruct.baseline)./exStruct.baseline;
exStruct.dF = (exStruct.baslineBleaching-exStruct.baseline)./exStruct.baseline;

end



