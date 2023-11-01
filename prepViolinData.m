function [data] = prepViolinData(dataCells, labels)

fullData = [];
dataLabels = [];
for i = 1:length(dataCells)
    try
        fullData = [fullData dataCells{i}];
    catch
        fullData = vertcat(fullData, dataCells{i});
    end
    dataLabels = [dataLabels repmat(labels(i),1, length(dataCells{i}))];
end

data.fullData = fullData;
data.dataLabels = dataLabels;
end