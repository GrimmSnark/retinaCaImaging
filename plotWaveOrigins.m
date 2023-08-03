function plotWaveOrigins(exStructPath)

%% load in exStruct
exStruct = load(exStructPath);
disp('Loaded in exStruct.mat')

try
    exStruct = exStruct.exStruct;
catch
    exStruct = exStruct.metaData;
end

%% load in images
BV_imagePath = dir([exStructPath(1:end-13) '*BV*tif']);
BV_image = read_Tiffs(fullfile(BV_imagePath.folder,BV_imagePath.name));

SD_imagePath = dir([exStructPath(1:end-13) '_dF_SD.tif']);
SD_image = read_Tiffs(fullfile(SD_imagePath.folder,SD_imagePath.name));

% downsample
BV_imageDS = imresize(BV_image, 0.25);

%% load into FIJI
stackImpBV = MIJ.createImage('BV', BV_imageDS, 1);
stackImpSD = MIJ.createImage('SD', SD_image, 1);

ij.process.ImageConverter(stackImpBV).convertToGray8();
ij.process.ImageConverter(stackImpSD).convertToGray8();

%% adjust brightness contrast and return to matlab
MIJ.run("Brightness/Contrast...");

% Sets up diolg box to allow for user input to choose cell ROIs
opts.Default = 'Continue';
opts.Interpreter = 'tex';

questText = [{'Manually adjust brightness/contrast'} ...
    {'If you are happy to move on with analysis click  \bfContinue\rm'} ...
    {'Or click  \bfExit Script\rm or X out of this window to exit script'}];

response = questdlg(questText, ...
    'Adjust contrast', ...
    'Continue', ...
    'Exit Script', ...
    opts);

% deals with repsonse
switch response

    case 'Continue' % if continue, goes on with analysis

    case 'Exit Script'
        return
end

MIJ.run("Merge Channels...", "c1=SD c2=SD c3=SD c6=BV create keep");
MIJ.run("RGB Color");


comImage = uint8(MIJ.getCurrentImage);

%% clean up
nImages = ij.WindowManager.getImageCount;
for i = 1: nImages
    stackImpWaves = ij.IJ.getImage();
    stackImpWaves.changes = false;
    stackImpWaves.close;
end

%% add wave start locations to image
figH = figure;
imshow(comImage, 'Border','tight');
hold on

for x = 1:length(exStruct.waves.centerPerFrame)
    waveStarts(x,:) = exStruct.waves.centerPerFrame{x}(1,:);
end

scatter(waveStarts(:,1),waveStarts(:,2),15,'g','filled');

%% save image
saveas(figH, [exStructPath(1:end-13) '_waveStarts.tif']);
close();
end