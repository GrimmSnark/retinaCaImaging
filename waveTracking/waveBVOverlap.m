function exStruct = waveBVOverlap(exStruct)

waveTable = exStruct.waves.waveTable;

for w = 1:max(waveTable.waveNumber)

    % get all objects for a single waves
    subTabIndx = waveTable.waveNumber == w;
    subTable = waveTable(subTabIndx,:);

    % get all the pixels involved in wave
    wavePixels = unique(cat(1,subTable.PixelIdxList{:}));
    waveIm = zeros(size(exStruct.waves.retinaBoundMask));
    waveIm(wavePixels) = 1;
    waveBound = bwboundaries(waveIm);
    wavePoly = polyshape(waveBound{1}(:,2),waveBound{1}(:,1));


    % BV poly and intersection
    try
        BV_poly = polyshape(exStruct.waves.BV_poly(:,1), exStruct.waves.BV_poly(:,2));
    catch
        exStruct.waves.BV_poly = exStruct.waves.BV_Position;
        exStruct.waves = rmfield(exStruct.waves, "BV_Position");
        BV_poly = polyshape(exStruct.waves.BV_poly(:,1), exStruct.waves.BV_poly(:,2));
    end
    waveBVIntersect = intersect(wavePoly,BV_poly);

    if waveBVIntersect.NumRegions > 0
        areaWave_pixel(w,1) = area(wavePoly);
        areaUnderBV_pixel(w,1) = area(waveBVIntersect);

%     figure
%     retBound = bwboundaries(exStruct.waves.retinaBoundMask');
%     retPoly = polyshape(retBound{1}(:,1),retBound{1}(:,2));
%     imshow(waveIm);
%     hold on
%     plot(retPoly);
%     set(gca,'YDir','reverse')
%     ylim([0 512]);
%     xlim([0 512]);
%     
%     hold on
%     plot(BV_poly);
%     plot(wavePoly);
%     plot(polyout);
%     close;

    else
        areaWave_pixel(w,1) = area(wavePoly);
        areaUnderBV_pixel(w,1) = 0;

    end
end

waveBVOverlap = table(areaWave_pixel, areaUnderBV_pixel );

waveBVOverlap.areaWave_micron = waveBVOverlap.areaWave_pixel * exStruct.downsampledRes;
waveBVOverlap.areaUnderBV_micron = waveBVOverlap.areaUnderBV_pixel * exStruct.downsampledRes;
waveBVOverlap.areaUnderBV_percent = (waveBVOverlap.areaUnderBV_pixel ./ waveBVOverlap.areaWave_pixel) * 100;


exStruct.waves.waveBVOverlap = waveBVOverlap;

end