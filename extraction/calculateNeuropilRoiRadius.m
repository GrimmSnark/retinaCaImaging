function radiusForNeuropil = calculateNeuropilRoiRadius(rois)
% Caclulates the Neuropil ROI radius based on the average width and height
% of the cell ROIs
%
% Input: cellular rois passed into MATLAB with rm.getRoisAsArray
%
% Output: radius in pixels for neuropil ROIs

%% Pre-Process ROIs and explicitly convert them to a ShapeRoi. 
% Sometimes ImageJ will use a different ROI type, and this code will throw
% an error as it expects the ROIs to inherit specific methods from the 
% ShapeROI class.

% Setup ROI Manager with MIJ
RC = ij.plugin.frame.RoiManager();
RM = RC.getInstance();

cellRois = rois;
for i =1:length(cellRois)
    cellName = cellRois(i).getName;
    cellRois(i) = ij.gui.ShapeRoi(cellRois(i));
    cellRois(i).setName(cellName);
end

%% Get bounding boxes for all ROIs

% Gets bounding info in rectangle object
for i=1:length(cellRois)
   bounding{i} = cellRois(i).getBounds; 
end

% Exacts to structure
for i =1:length(bounding)
   boundingStruct(i) = struct(bounding{i}); 
end

%% Gets averages

widthAverage = mean([boundingStruct(:).width]);
hightAverage = mean([boundingStruct(:).height]);

radiusForNeuropil = round(mean([widthAverage hightAverage]))/2;

end
