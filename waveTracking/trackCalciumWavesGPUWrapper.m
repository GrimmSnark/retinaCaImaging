function trackCalciumWavesGPUWrapper(dataDir, fileStartNo)
% This function runs trackCalciumWaves on a full folder of recordings which
% have the dF_F movies created on the GPU. This is the preferred version...
%
% See also: trackCalciumWavesGPU
%
% Input- folderPath- folder path which contains .nd2 files from FLAME 
%                    system
%


if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

dFPaths = dir([dataDir '*dF_F.tif']);

for i = fileStartNo:length(dFPaths)
    disp(['On experiment No. ' num2str(i) ' of '  num2str(length(dFPaths))]);
    trackCalciumWavesGPU(fullfile(dFPaths(i).folder, dFPaths(i).name));
end

end