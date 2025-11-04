function prepRetinaCalciumWrapper(folderPath, fileStartNo, motionCorrFlag,  motionCorrectionType,  createDFPixelMovieFlag, zeroDFStack,  channelOrg)
% This function runs prepRetinaCalcium on a full folder of recordings
% See also: prepRetinaCalcium
%
% Input- folderPath- folder path which contains .nd2 files from FLAME 
%                    system
%
%        fileStartNo - File number to start the wrapper on, DEFAULT == 1
%
%        motionCorrFlag- Flag do motion correction on raw images
%
%        motionCorrectionType - DFT-based subpixel method
%                               ('subMicronMethod')
%                             - non-rigid NoRM Corr registration
%                               ('nonRigid')
%
%        createDFPixelMovieFlag: flag for creating pixel wise DF_F stack
%        and takes up time/space 0 = not saved, 1 = saved (DEFAULT)
%
%        zeroDFStack - 0/1 to make and negative DF signal zero for the DF
%                      movies DEFAULT = 1
%
%        channelOrg - two element vector showing organisation of calcium 
%                     then blood vessel channels, first element is the
%                     calcium channel no, the second is the blood vessel
%                     channel DEFAULT [1 2]

%% defaults

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

if nargin < 3 || isempty(motionCorrFlag)
    motionCorrFlag = 1;
end

if nargin < 4 || isempty(motionCorrectionType)
    motionCorrectionType = 'subMicronMethod';
end

if nargin < 5 || isempty(createDFPixelMovieFlag)
    createDFPixelMovieFlag = 1;
end

if nargin < 6 || isempty(zeroDFStack)
    zeroDFStack = 1;
end


if nargin <7 || isempty(channelOrg)
    channelOrg = [1 2];
end

files = dir([folderPath '**\*rec*.nd2']);

% remove snaps
files = files(~contains({files(:).name},'snap'));

if isempty(files)
    files = dir([folderPath '*.czi']);
end

if isempty(files)
    files = dir([folderPath '*.tif']);
end

for i = fileStartNo:length(files)
    disp(['Preprocessing file no. ' num2str(i) ' of ' num2str(length(files))])
    prepRetinaCalcium(fullfile(files(i).folder,files(i).name),motionCorrFlag,motionCorrectionType,createDFPixelMovieFlag,zeroDFStack, channelOrg);
end

end