function plotRelativeDisWaveStarts_Wrapper(exStructPath, fileStartNo, useDiologue)

if nargin < 1 || isempty(exStructPath)
    exStructPath= uigetdir([],...
        'Image File Selector');
end

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

if nargin < 3 || isempty(useDiologue)
    useDiologue = 1;
end


files = dir([exStructPath '\*.mat']);

for i = fileStartNo:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
     plotRelativeDisWaveStarts(fullfile(files(i).folder,files(i).name), useDiologue);
end

end