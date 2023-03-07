function plotCellDF(exStructPath)
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


exStruct = load(exStructPath);
exStruct = exStruct.exStruct;
cells = exStruct.cells;

%% create folder

saveDir = [exStruct.filePath(1:end-4) '_cellPlots'];

if ~exist(saveDir,"dir")
    mkdir(saveDir);
end

%% run through the cells to plot

timeBase = (1:size(cells.dF,2))/exStruct.rate; % time base in seconds

for c = 1:cells.cellCount
figH = figure('units','normalized','outerposition',[0 0 1 1]);

plot(timeBase, cells.dF(c,:));

title(['DF for Cell No' num2str(c)]);

xlabel('Time in Sec');
xlim([0 timeBase(end)]);

ylabel('DF');
tightfig();

saveas(figH,fullfile(saveDir, ['CellDF_' num2str(c) '.png']));
close

end
end