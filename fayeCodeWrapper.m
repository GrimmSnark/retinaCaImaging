function fayeCodeWrapper(path)

exStructPaths = dir([path '\**\*exStruct.mat']);

disp(['Found ' num2str(length(exStructPaths)) ' experiment folders, moving on....'])

for i =1:length(exStructPaths)
    currentExPath = fullfile(exStructPaths(i).folder, exStructPaths(i).name);

    fayeCode(currentExPath);

    % uncomment this if you want to close every figure after each mat file processed
    % close
end
end