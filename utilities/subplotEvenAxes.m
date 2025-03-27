function subplotEvenAxes(figH, axis2Use, subs2Align)
% Sets both X,Y,Z axes to be the same across an entire subplot image
%
% Inputs: figH - figure handle to even out
%         axis2Use - vector to specify axes to even out (x y z)
%                    DEFAULT = all axis [1 1 1]

%% defaults
if nargin < 2 || isempty(axis2Use)
    axis2Use = [1 1 1];
end

%% check if the handle is a figure
if isgraphics(figH, 'axes')

    axH = figH;

elseif ishandle(figH) && findobj(figH,'type','figure')==figH

    axH2 = findall(figH,'type','axes');
end

if nargin < 3 || isempty(subs2Align)
    subplotInd = ones(size(axH));
else
    subplotInd = false(size(axH));
    subplotInd(subs2Align) = 1;
end

axH = axH(subplotInd);

%% X axis
if axis2Use(1) == 1
    axXLim =  get(axH,'XLim');
    maxX = max([axXLim{:}]);
    minX = min([axXLim{:}]);
    set(axH,'XLim',[minX maxX]);
end

%% Y axis
if axis2Use(2) == 1
    axYLim =  get(axH,'YLim');
    maxY = max([axYLim{:}]);
    minY = min([axYLim{:}]);
    set(axH,'YLim',[minY maxY]);
end

%% Z axis
if axis2Use(3) == 1
    try
        axZLim =  get(axH,'ZLim');
        maxZ = max([axZLim{:}]);
        minZ = min([axZLim{:}]);
        set(axH,'ZLim',[minZ maxZ]);
    catch
        
    end
end
end