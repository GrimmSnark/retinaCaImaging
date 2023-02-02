function prepRetinaCalciumWrapper(folderPath)

files = dir([folderPath '*.nd2']);
for i = 1:length(files)
    disp(['On file no. ' num2str(i) ' of ' num2str(length(files))])
    prepRetinaCalcium(fullfile(files(i).folder,files(i).name));
end

end