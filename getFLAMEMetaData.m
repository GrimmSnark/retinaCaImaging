function metaData = getFLAMEMetaData(xmlStruct)
% Grabs the metadata from the FLAME system into matlab format


%% get all the metadata

% get objective info
objectiveMag = double(xmlStruct.getObjectiveCalibratedMagnification(0,0));
objectiveNA = double(xmlStruct.getObjectiveLensNA(0,0));
objectiveModel = string(xmlStruct.getObjectiveModel(0,0));

% get imageID
imageName = string(xmlStruct.getImageName(0));

% get resolution
imageSizePixels = xmlStruct.getPixelsSizeX(0).getValue;
pixelSize = double(xmlStruct.getPixelsPhysicalSizeX(0).value);
pixelSizeUnit = string(xmlStruct.getPixelsPhysicalSizeX(0).unit.getSymbol);

% get channel num
chanNum = xmlStruct.getChannelCount(0);

% get colour channel info
for i =1:chanNum
    emWave(i) = double(xmlStruct.getChannelEmissionWavelength(0,i-1).value);
    emWaveName{i} = string(xmlStruct.getChannelName(0,i-1));
end

imNo = xmlStruct.getPlaneCount(0);

% get frame times
for i = 1:imNo
    frameTime(i) = double(xmlStruct.getPlaneDeltaT(0,i-1).value);
    exposureTime(i) = double(xmlStruct.getPlaneExposureTime(0,i-1).value);
    imgColIndx(i) = xmlStruct.getPlaneTheC(0,i-1).getValue;
    imgTimepointIndx(i) = xmlStruct.getPlaneTheT(0,i-1).getValue;
    imgZPos(i) = double(xmlStruct.getPlanePositionZ(0,i-1).value);
    imgZIndx(i) = xmlStruct.getPlaneTheZ(0,i-1).getValue;
end

framePeriod = mean(diff(unique(frameTime)));
fps = 1/framePeriod;

%% compile into struct

metaData.imageName = imageName;

% objective info
metaData.objective.objectiveMag = objectiveMag;
metaData.objective.objectiveNA = objectiveNA;
metaData.objective.objectiveModel = objectiveModel;

% image specs
metaData.image.pixelNum = imageSizePixels;
metaData.image.pixelSize = pixelSize;
metaData.image.pixelSizeUnit = pixelSizeUnit;

% image colors
metaData.colours.emWavelength = emWave;
metaData.colours.emWavelengthName = emWaveName;

% frame info
metaData.frameInfo.frameTime = frameTime;
metaData.frameInfo.exposureTime = exposureTime;
metaData.frameInfo.imgColourIndx = imgColIndx;
metaData.frameInfo.imgTimepointIndx = imgTimepointIndx;
metaData.frameInfo.imgZPos = imgZPos;
metaData.frameInfo.imgZIndx = imgZIndx;

metaData.framePeriod = framePeriod;
metaData.fps = fps;

end