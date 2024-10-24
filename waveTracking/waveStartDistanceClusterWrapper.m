function waveStartDistanceClusterWrapper(filePath)

exStructPaths = dir([filePath '\**\*.mat']);


for i = 1:length(exStructPaths)
    if ~contains(exStructPaths(i).folder,'bad')
        load(fullfile(exStructPaths(i).folder, exStructPaths(i).name));

        try
            disp(['On file ' num2str(i) ' of ' num2str(length(exStructPaths))]);
            exStruct = waveStartDistanceCluster(exStruct);
            save(fullfile(exStructPaths(i).folder, exStructPaths(i).name), "exStruct", '-v7.3');

        catch

        end
    end

end

end