function BV_D1_2 = calculateVascularD1_2(exStructPath)
% creates the blood vessel D1/2 metric for each recording by averaging the
% D1/2 for all blood vessel vertices 
% Already implemented in plotRelativeDisWaveStarts for ongoing experiments.
% This is only for first run through/plotting for ES. 

%% load in exStruct
load(exStructPath);
disp('Loaded in exStruct.mat')

%% load in polygons

% optic nerve poly
centroidStruct = regionprops(exStruct.waves.opticNerveMask,"Centroid");
opticNerveCentroid = centroidStruct.Centroid;

% retina boundary poly
retinaBoundPoly = bwboundaries(exStruct.waves.retinaBoundMask');
retinaBoundShape = polyshape(retinaBoundPoly{1});

% blood vessel boundary with added error checking FFS
try
    BV_poly = exStruct.waves.BV_poly;
catch
    exStruct.waves.BV_poly = exStruct.waves.BV_Position;
    exStruct.waves = rmfield(exStruct.waves, 'BV_Position');
    BV_poly = exStruct.waves.BV_poly;
end

if isempty(exStruct.waves.BV_poly)
    exStruct.waves.BV_poly = exStruct.waves.BV_Position;
    exStruct.waves = rmfield(exStruct.waves, 'BV_Position');
    BV_poly = exStruct.waves.BV_poly;
end

%% Run through BV_poly to get D1/2s

for i=1:length(BV_poly)

    % fit line through retinal centre and BV vertex position through to
    % edge of image
    pFit(i,:,:) = fitStraightLine(opticNerveCentroid, BV_poly(i,:), [0 512]);

    % get both side intersections on the line
    currLine = squeeze(pFit(i,:,:));
    [inR, outR] = intersect(retinaBoundShape,  currLine);

    % find the shortest path between BV and retina bound
    [~,indUsed] = pdist2(inR, BV_poly(i,:),'euclidean','Smallest',1);
    retinaCrossingCoor = inR(indUsed,:);
    line2Use = [opticNerveCentroid; retinaCrossingCoor];

    center2RetinaEdge = pdist2(retinaCrossingCoor, opticNerveCentroid,'euclidean','Smallest',1);
    center2BV_Vert = pdist2(BV_poly(i,:), opticNerveCentroid,'euclidean','Smallest',1);

    vertexD1_2(i) = center2BV_Vert/center2RetinaEdge;
end

BV_D1_2 = mean(vertexD1_2);
exStruct.waves.BV_D1_2 = BV_D1_2;

save(exStructPath, "exStruct", '-v7.3');

end

