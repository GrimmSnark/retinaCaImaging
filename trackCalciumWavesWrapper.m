function trackCalciumWavesWrapper(folderPath)

files = dir([folderPath '*dF_F.tif']);

for i = 1:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
    try
        trackCalciumWaves(fullfile(files(i).folder,files(i).name), 0.98);
    catch
    end
end

end