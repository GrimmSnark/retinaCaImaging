function plotWaveOrigins_Wrapper(exStructPath, fileStartNo)

if nargin < 1 || isempty(exStructPath)
    exStructPath= uigetdir([],...
        'Image File Selector');
end

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

files = dir([exStructPath '\*.mat']);

for i = fileStartNo:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
     plotWaveOrigins(fullfile(files(i).folder,files(i).name));
end

end