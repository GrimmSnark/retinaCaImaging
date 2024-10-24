function errorStack = fayeCodeWrapper(path)

if nargin < 1 || isempty(path)
    [path] = uigetdir([],...
        'Folder to Process Selector');
end

exStructPaths = dir([path '\**\*exStruct.mat']);

disp(['Found ' num2str(length(exStructPaths)) ' experiment folders, moving on....'])
count = 1;
for i =1:length(exStructPaths)
    currentExPath = fullfile(exStructPaths(i).folder, exStructPaths(i).name);

    try
    fayeCode(currentExPath);

    % uncomment this if you want to close every figure after each mat file processed
    % close

    catch ME
        text2Display = ['Unable to process ' currentExPath ' at ' ME.stack(1).file ' line ' num2str(ME.stack(1).line) ];
        warning(text2Display);

        errorStack{count} = text2Display;
        count = count +1;

    end
end
end