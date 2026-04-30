function [imStack, imageMetaData] = readFLAMEData(filePath)
% read in the .nd2 file from the FLAME system and make metadata structure
%
% Inputs: filePath - fullfile to nd2 file
%
% Outputs: imStack - image stack from the nd2 file
%
%          metaData - meta data structure extracted from nd2 file

%% read in nd2 file
try
    tic
    imageStruct = bfopen2_parallel(filePath);
    toc
catch
    disp('Deleting bfmemo file and retrying');
    [fparts,name,ext] = fileparts(filePath);
    bfmemoFilepath = dir(fullfile(fparts, ['.' name ext '.bfmemo']));
    delete(fullfile(bfmemoFilepath.folder, bfmemoFilepath.name));

    imageStruct = bfopen2_parallel(filePath);
    toc
end

%% get the meta data
% try
    omeMeta = imageStruct{1, 4};
    hashMeta = imageStruct{1, 2};
    imageMetaData = getFLAMEMetaData(omeMeta, hashMeta);
    imageMetaData.filePath = filePath;
% catch
%     clear javaclasspath
%     [imageMetaData] = bfinfo(filePath);
%     imageMetaData = metaData{1};
%     imageMetaData.filePath = filePath;
% end
%% get the images
imStack = imageStruct{1,1}(:,1);

% split into channels
imStack = cat(3,imStack{:});
end