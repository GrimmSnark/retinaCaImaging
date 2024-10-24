function CaAnalysisFLAME_Wrapper(folderPath, fileStartNo, baselineSubType)
% This function runs calcium cell analysis on a full folder of recordings
% See also: CaAnalysisFLAME
%
% Input- folderPath- folder path which contains .nd2 files from FLAME 
%                    system
%
%        fileStartNo - File number to start the wrapper on, DEFAULT == 1
%
%        baselineSubType - Switch case for calcium imaging baseline
%                          subtraction type DEFAULT == 2
%                          1: Expotential bleaching fit then rolling ball
%                             baseline median percentile filter (use for
%                             highly packed cells, ie retina)
%                          2: Annulus neuropil subtraction and kernal 
%                             density estimation for percentile filter (
%                             use for loosely packed cells, ie culture or
%                             in-vivo brain)
%% Defaults
if nargin < 1 || isempty(folderPath)
    folderPath= uigetdir([],...
        'Image Folder Selector');
end

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

if nargin < 3 || isempty(baselineSubType)
    baselineSubType = 2;
end

%%
files = dir([folderPath '\*.mat']);

for i = fileStartNo:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
    CaAnalysisFLAME(fullfile(files(i).folder,files(i).name),baselineSubType);
end

end