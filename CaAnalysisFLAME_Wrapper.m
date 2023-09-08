function CaAnalysisFLAME_Wrapper(exStructPath, baselineSubType)

if nargin < 1 || isempty(exStructPath)
    exStructPath= uigetdir([],...
        'Image File Selector');
end

if nargin < 2 || isempty(baselineSubType)
    baselineSubType = 2;
end

files = dir([exStructPath '*.mat']);

for i = fileStartNo:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
    CaAnalysisFLAME(fullfile(files(i).folder,files(i).name),baselineSubType);
end

end