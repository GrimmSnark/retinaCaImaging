function parentDir =  getParentDirFull(filePath)

parts = strsplit(filePath, '/');

if length(parts)<2
    parts = strsplit(filePath, '\');
end

parentDir = strjoin(parts(1:end-1), filesep);

end