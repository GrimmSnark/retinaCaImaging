function trackCalciumWavesGPUWrapper(dataDir, fileStartNo)

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

dFPaths = dir([dataDir '*dF_F.tif']);

for i = fileStartNo:length(dFPaths)
    trackCalciumWavesGPU(fullfile(dFPaths(i).folder, dFPaths(i).name),1);
end

end