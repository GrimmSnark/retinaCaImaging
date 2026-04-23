function calculateDrugResponseAUC(exStructPath, baselineStartTime, drugResponseStartTime, chunkSizeTime)
% Calculates Area Under Curve for a period of baseline and drug response
% based on inputs. Adds info into exStruct and creates an excel file
%
% Inputs: exStructPath - filepath for exStruct.mat to process
%         baselineStartTime - time in seconds for start of baseline AUC
%         drugResponseStartTime - time in seconds for start of drug AUC
%         chunkSizeTime - AUC chunk size in sec, DEFAULT == 60s

%% load exStruct

if nargin < 1 || isempty(exStructPath)
    [file, path] = uigetfile({'*.mat'},...
        'Image File Selector');

    exStructPath = fullfile(path,file);
end

if nargin < 4 || isempty(chunkSizeTime)
    chunkSizeTime = 60; % seconds
end

exStruct = load(exStructPath);
exStruct = exStruct.exStruct;
fps = round(exStruct.fps);

chunkSizeFrames = round(chunkSizeTime * fps);

%% get baseline AUC
baselineStartFrames = baselineStartTime * fps;

for c= 1:exStruct.cellCount
    baselineDFwindow = exStruct.cells.dF(c, baselineStartFrames:baselineStartFrames+chunkSizeFrames);

    % get everything above zero
    yData = baselineDFwindow(baselineDFwindow>0);
    xData = find(baselineDFwindow>0);
    % Calculate the area under the dFTrace for each spike
    if length(yData) > 1
        AUC_baseline(c) = trapz(xData,yData);
    else
        AUC_baseline(c) = 0;
    end
end

AUC_baseline =AUC_baseline';

%% get drug AUC
drugResponseStartFrames = drugResponseStartTime * fps;

for c= 1:exStruct.cellCount
    baselineDFwindow = exStruct.cells.dF(c, drugResponseStartFrames:drugResponseStartFrames+chunkSizeFrames);

    % get everything above zero
    yData = baselineDFwindow(baselineDFwindow>0);
    xData = find(baselineDFwindow>0);
    % Calculate the area under the dFTrace for each spike
    if length(yData) > 1
        AUC_drug(c) = trapz(xData,yData);
    else
        AUC_drug(c) = 0;
    end
end

AUC_drug =AUC_drug';

%% add into exStruct and save
exStruct.cells.AUC_drugResponse.baselineStartT = baselineStartTime;
exStruct.cells.AUC_drugResponse.drugStartT = drugResponseStartTime;
exStruct.cells.AUC_drugResponse.windowLen = chunkSizeTime;

exStruct.cells.AUC_drugResponse.baseline = AUC_baseline;
exStruct.cells.AUC_drugResponse.drug = AUC_drug;

save(exStructPath, "exStruct", '-v7.3');

%% save the output as excel

AUC_table = table(AUC_baseline, AUC_drug );

writetable(AUC_table, [exStruct.filePath(1:end-4) '_AUC_drug.xlsx']);


end