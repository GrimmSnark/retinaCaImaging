function prepRetinaCalciumWrapper(folderPath, fileStartNo)

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

files = dir([folderPath '*.nd2']);

if isempty(files)
    files = dir([folderPath '*.czi']);
end

for i = fileStartNo:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
    prepRetinaCalcium(fullfile(files(i).folder,files(i).name),[],[],[],1);
end

end