function plotWaveProbabilityIm(exStructPath)

%% load in exStruct
exStruct = load(exStructPath);
exStruct = exStruct.exStruct;
disp('Loaded in exStruct.mat')

%% load in waveCol movie

[saveFolder, name] = fileparts(exStructPath);
waveColPath = fullfile(saveFolder, [name(1:end-9) '_waveCol.tif']);
waveIm = read_Tiffs(waveColPath,[],1);

%% 
% remove white from RGB stacks

for i =1:size(waveIm,4)

    tempIm = waveIm(:,:,:,i);

    % seperate into RGB
    tempImR = tempIm(:,:,1);
    tempImG = tempIm(:,:,2);
    tempImB = tempIm(:,:,3);

    % find the indexes per channel
    whitePixR = find(tempImR==255);
    whitePixG = find(tempImG==255);
    whitePixB = find(tempImB==255);

    % group together and see how many are in all three channels
    allMaxPixels = [whitePixR;whitePixG;whitePixB];
    [counts,gValues] = groupcounts(allMaxPixels);

    whiteIndexs = gValues(counts==3);

    % convert to row, col values
    [r, c] = ind2sub(size(tempImR), whiteIndexs);


    % remove white pixels
    for x =1:length(r)
        tempIm(r(x),c(x),:) = [0 0 0];
    end

    waveIm(:,:,:,i) = tempIm;
    binaryImStack(:,:,i) = im2bw(tempIm);
end


probWaveIm = (sum(binaryImStack,3))/size(binaryImStack,3);

% rescaled
probWaveImRescale = rescale(probWaveIm, 0, 255);
% Typecasted to uint8 
probWaveImRescale = uint8(probWaveImRescale);
figH = imshow(probWaveImRescale,lcs,Border="tight");
hold on

% plot retina bound
retinaBounds = bwboundaries(exStruct.waves.retinaBoundMask);
plot(retinaBounds{1}(:,2),retinaBounds{1}(:,1),'Color','k', LineWidth=1);

% plot blood vessels
try
    BVPoints = [exStruct.waves.BV_poly; exStruct.waves.BV_poly(1,:)] ;
    plot(BVPoints(:,1), BVPoints(:,2),'Color','g','LineWidth',1);

catch
    BVPoints = [exStruct.waves.BV_Position; exStruct.waves.BV_Position(1,:)] ;
    plot(BVPoints(:,1), BVPoints(:,2),'Color','g','LineWidth',1);

end

% save
saveas(figH, fullfile(saveFolder, [name(1:end-9) '_waveProbMap.tif']));

close;
end