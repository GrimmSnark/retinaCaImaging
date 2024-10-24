function fixWaveAreaWrapper(folder2Search, startNo)

if nargin<2 || isempty(startNo)
    startNo = 1;
end


exStructPaths = dir([folder2Search '\**\*\*exStruct.mat']);

for i = startNo:length(exStructPaths)

    disp(['Fixing file ' num2str(i) ' of ' num2str(length(exStructPaths))])
    load(fullfile(exStructPaths(i).folder, exStructPaths(i).name));

    if isfield(exStruct, 'waves')

        exStruct = fixWaveArea(exStruct);
        save(fullfile(exStructPaths(i).folder, exStructPaths(i).name), "exStruct", '-v7.3');
    end
end

end