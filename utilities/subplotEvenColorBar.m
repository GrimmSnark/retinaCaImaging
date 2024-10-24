function subplotEvenColorBar(axH)
% Sets Colorbar limits to be the same across an entire subplot image

%% check if the handle is a figure
if ishandle(axH) && findobj(axH,'type','figure')==axH
    axH = findall(axH,'type','axes');
end

%% Set Color axis limits
axCLim =  get(axH,'Clim');
maxC = max([axCLim{:}]);
minC = min([axCLim{:}]);
set(axH,'CLim',[minC maxC]);

end