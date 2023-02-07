function createDFPixelMovieZeroed(exStructPath)
%% load in exStruct
exStruct = load(exStructPath);
disp('Loaded in exStruct.mat')

try
    exStruct = exStruct.exStruct;
catch
    exStruct = exStruct.metaData;
end

%% create dF Zero movie
dFZero = exStruct.dF;
dFZero(dFZero<0) = 0;

dFZero = uint16(mat2gray(dFZero) * 65535);

saveastiff(dFZero, [exStructPath(1:end-13) '_dF_F.tif']);
end