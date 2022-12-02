function [imStack, metaData] = readFLAMEData(filePath)
% read in the .nd2 file from the FLAME system and make metadata structure
%
% Inputs: filePath - fullfile to nd2 file
%
% Outputs: imStack - image stack from the nd2 file
%
%          metaData - meta data structure extracted from nd2 file

%% read in nd2 file
imageStruct = bfopen2(filePath);

%% get the meta data
omeMeta = imageStruct{1, 4};
metaData = getFLAMEMetaData(omeMeta);
metaData.filePath = filePath;

%% get the images
imStack = imageStruct{1,1}(:,1);

% split into channels
imStack = cat(3,imStack{:});
end