function plotWaveProbabilityImWrapper(folderPath, fileStartNo)

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

try
    files = dir([folderPath '**\*.mat']);
catch
    files = dir([folderPath '\**\*.mat']);
end

% remove snaps
files = files(~contains({files(:).folder},'bad'));

for i = fileStartNo:length(files)
    disp(['Plotting experiment no. ' num2str(i) ' of ' num2str(length(files))])
    plotWaveProbabilityIm(fullfile(files(i).folder,files(i).name));

end


end