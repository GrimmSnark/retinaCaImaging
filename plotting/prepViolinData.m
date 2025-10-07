function [data] = prepViolinData(dataCells, labels)
%% Prepares data for use in box plots and violin plots
% 
% Inputs: 
%           dataCells - Cell array of metrics or variables to unpack
%
%           labels - String array of length dataCells x label 
% 
% Output:   data- structure containing
%                      data.fullData- unrolled data into a vector
%
%                      data.dataLabels - table of data labels for each
%                                        vector entry in fullData
%
%                      data.uniqueLabelCombs - table of unique label
%                                              combinations
%
%                      data.uniqueLabelIndx - indexes for each data entry 
%                                             for the unique combinations 


%% Get the data unzipped

fullData = [];
dataLabels = table();
for i = 1:length(dataCells)
    try
        fullData = [fullData dataCells{i}];
    catch
        fullData = vertcat(fullData, dataCells{i});
    end

        dataLabels = [dataLabels; cell2table(cellstr(repmat(labels(i,:), length(dataCells{i}),1)))];
end

[uniqueLabelCombs,~, uniqueLabelIndx] = unique(dataLabels,"rows");

%% Pull together data
data.fullData = fullData;
data.dataLabels = dataLabels;
data.uniqueLabelCombs = uniqueLabelCombs;
data.uniqueLabelIndx = uniqueLabelIndx;
end