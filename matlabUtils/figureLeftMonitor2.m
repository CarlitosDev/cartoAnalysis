function f = figureLeftMonitor2()

  % Set the figure 1st monitor on the left
  f = figure();

% reposition to left half of the second monitor
monitorPositions = get(0,'MonitorPositions');
[numMonitors, ~] = size(monitorPositions);

if numMonitors > 1
  
  screenSize       = monitorPositions(2,:);
  screenWidth      = screenSize(3)-screenSize(1)+1;
  screenHeight     = screenSize(4);


  % to-do...make sure the 2nd one is smaller than the 1st one
  secondMonitorOffset = monitorPositions(1,4)-monitorPositions(2,4);

  currentFigure  = gcf();
  figurePosition = get(currentFigure, 'OuterPosition');



  fromLeft   = screenSize(1);
  fromBottom = 1+secondMonitorOffset;
  newWidth   = screenWidth/2;
  newHeight  = screenHeight;
  set(f, 'units', 'pixels', 'OuterPosition',  [fromLeft fromBottom newWidth newHeight]);
  
end