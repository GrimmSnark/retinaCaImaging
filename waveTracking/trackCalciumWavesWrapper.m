function trackCalciumWavesWrapper(folderPath)
% This function runs trackCalciumWaves on a full folder of recordings which
% have the dF_F movies created on the CPU, can be slow
%
% See also: trackCalciumWaves
%
% Input- folderPath- folder path which contains .nd2 files from FLAME 
%                    system
%

files = dir([folderPath '*dF_F.tif']);

for i = 1:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
    try
        trackCalciumWaves(fullfile(files(i).folder,files(i).name));
    catch
    end
end

end