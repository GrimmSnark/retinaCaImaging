function loadImagingData2FIJI(filepath2Use, registerFlag)
% loads processed data into ImageJ and preforms motion correction if you
% want
%
% Inputs: filepath2Use - folder for processed data to load
%
%         registerFlag - 0/1 to preform motion correction
%                        DEFAULT = 1, do correction


%% set defaults

if nargin <2 || isempty(registerFlag)
   registerFlag = 1; 
end

usingRawData = 0;
intializeMIJ;

%% gather data
% check whether filepath is processed folder
folder2Process = dir([filepath2Use '\**\experimentStructure.mat']);

% tries to fix if we happened to use the raw data path
if isempty(folder2Process)
    filePath = createSavePath(filepath2Use,1,1);
    folder2Process = dir([filePath '\**\experimentStructure.mat']);
end

if isempty(folder2Process)
    folder2Process = filepath2Use;
    usingRawData  =1;
end

if usingRawData ==0
    
    % load experimentStructure
    load([folder2Process.folder '\experimentStructure.mat']);
    
    % read in tiff file
    vol = read_Tiffs(experimentStructure.fullfile,1);
    if ndims(vol) ~=3
        vol = readMultipageTifFiles(experimentStructure.prairiePath);
    end
    
    
    % check number of channels in imaging stack
    channelIndxStart = strfind(experimentStructure.filenamesFrame{1}, '_Ch');
    for i =1:length(experimentStructure.filenamesFrame)
        channelIdentity{i} = experimentStructure.filenamesFrame{i}(channelIndxStart:channelIndxStart+3);
    end
    channelNo = unique(channelIdentity);
    
    % chooses correct channel to analyse in multichannel recording
    if length(channelNo)>1
        volSplit =  reshape(vol,size(vol,1),size(vol,2),[], length(channelNo));
        vol = volSplit;
    end
    
    if registerFlag == 1
    for q = 1:size(vol,4)
        % apply imageregistration shifts if there are shifts to apply
        if isprop(experimentStructure, 'options_nonrigid') && ~isempty(experimentStructure.options_nonrigid) % if using non rigid correctionn
            registeredVol(:,:,:,q) = apply_shifts(vol(:,:,:,q),experimentStructure.xyShifts,experimentStructure.options_nonrigid);
        elseif  ~isempty(experimentStructure.xyShifts)
            registeredVol(:,:,:,q) = shiftImageStack(vol(:,:,:,q),experimentStructure.xyShifts([2 1],:)'); % Apply actual shifts to tif stack
        else % if there are no motion correction options, ie the image stack loaded is already motion corrected
            registeredVol(:,:,:,q) = vol(:,:,:,q);
        end
    end
    else
        registeredVol = vol;
    end
    
    if length(channelNo)>1
        registeredVol = reshape(registeredVol, size(registeredVol,1),size(registeredVol,2),[],1);
    end
else
    
end

% transfers to FIJI
registeredVolMIJI = MIJ.createImage( 'Registered Volume', registeredVol,true);




end