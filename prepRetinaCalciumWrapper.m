function prepRetinaCalciumWrapper(folderPath, fileStartNo, channelOrg)

if nargin < 2 || isempty(fileStartNo)
    fileStartNo = 1;
end

if nargin <3 || isempty(channelOrg)
    channelOrg = [];
end

files = dir([folderPath '**\*rec*.nd2']);

% remove snaps
files = files(~contains({files(:).name},'snap'));

if isempty(files)
    files = dir([folderPath '*.czi']);
end

for i = fileStartNo:length(files)
    disp(['Preprocessing file no. ' num2str(i) ' of ' num2str(length(files))])
    prepRetinaCalcium(fullfile(files(i).folder,files(i).name),[],[],[],1, channelOrg);
end

end