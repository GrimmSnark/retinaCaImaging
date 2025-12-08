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

channel2Use = 1;
usingRawData = 0;
intializeMIJ;

%% gather data
% check whether filepath is processed folder
folder2Process = dir([filepath2Use '\**\*_ExStruct.mat']);

% % tries to fix if we happened to use the raw data path
% if isempty(folder2Process)
%     filePath = createSavePath(filepath2Use,1,1);
%     folder2Process = dir([filePath '\**\experimentStructure.mat']);
% end
%
if isempty(folder2Process)
    folder2Process = filepath2Use;
    % usingRawData  =1;
end

if usingRawData ==0

    % load experimentStructure
    load([folder2Process]);

    if ~exist(exStruct.filePath)
        parentDir =  getParentDirFull(folder2Process);

        imPath = dir([parentDir '\*.nd2']);
        imPath = fullfile(imPath.folder,imPath.name);
        exStruct.filePath = imPath;
    end

    [~, ~, ext ] = fileparts(exStruct.filePath);

    if usingRawData ==0
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

    end

else

end

% transfers to FIJI
registeredVolMIJI = MIJ.createImage( 'Registered Volume', registeredVol,true);




end