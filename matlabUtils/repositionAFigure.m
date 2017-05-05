% reposition to left half of the main monitor
monitorPositions = get(0,'MonitorPositions');
screenSize       = get(0,'ScreenSize');
screenWidth      = screenSize(3);
screenHeight     = screenSize(4);



% Set the figure 1st monitor on the left
f = figure();
currentFigure  = gcf();
figurePosition = get(currentFigure, 'OuterPosition');

% Let's assume is windows 7
if ispc
  taskBarOffset = 40;
else
  taskBarOffset = 0;
end


fromLeft   = 1;
fromBottom = 1 + taskBarOffset;
newWidth   = screenWidth/2;
newHeight  = screenHeight - fromBottom;
set(f, 'units', 'pixels', 'OuterPosition',  [fromLeft fromBottom newWidth newHeight])


% Set the figure 1st monitor on the left
f2 = figure();
currentFigure  = gcf();
figurePosition = get(currentFigure, 'OuterPosition');

if ispc
  taskBarOffset = 40;
else
  taskBarOffset = 0;
end


fromLeft   = 1 + screenWidth/2;
fromBottom = 1 + taskBarOffset;
newWidth   = screenWidth/2;
newHeight  = screenHeight - fromBottom;
set(f2, 'units', 'pixels', 'OuterPosition',  [fromLeft fromBottom newWidth newHeight])


