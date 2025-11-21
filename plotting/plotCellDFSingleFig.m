function plotCellDFSingleFig(exStructPath, offset)
% This function runs a very basic plotting of all cell DFs calcluated
% during the calcium analysis pipeline
%
% Written by Michael Savage (michael.savage2@ncl.ac.uk)
%
% Inputs-  exStructPath: filepath for exStruct.mat to process, or leave
%                        empty to start a GUI file picker

%% load exStruct

if nargin < 1 || isempty(exStructPath)
    [file, path] = uigetfile({'*.mat'},...
        'Image File Selector');

    exStructPath = fullfile(path,file);
end

if nargin <2 || isempty(offset)

    offset = 0.4;
end


exStruct = load(exStructPath);
exStruct = exStruct.exStruct;
cells = exStruct.cells;

%% create folder

%% run through the cells to plot

try
    timeBase = (1:size(cells.dF,2))/exStruct.rate; % time base in seconds
catch
    timeBase = (1:size(cells.dF,2))/exStruct.fps; % time base in seconds
end

figH = figure('units','normalized','outerposition',[0 0 1 1]);
hold on

plotOffset = 0;
for c = 1:cells.cellCount

    plotOffset = plotOffset + offset;
    plot(timeBase, cells.dF(c,:) + plotOffset);
    ytickPos(c) = plotOffset;

    % title(['DF for Cell No' num2str(c)]);

    xlabel('Time in Sec');
    xlim([0 timeBase(end)]);

    ylabel('DF');
end

yticks(ytickPos);
yticklabels(1:c);
tightfig();
saveas(figH,[fullfile(exStruct.filePath(1:end-4)) '_CellDFs.png']);

end