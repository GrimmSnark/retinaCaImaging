function createCaWaveRasterPlot_Wrapper(filePath)

files = dir([filePath '\**\*ExStruct*']);

for i = 1:length(files)
    try
        createCaWaveRasterPlotV2(fullfile(files(i).folder, files(i).name));
    catch
    end
end

end