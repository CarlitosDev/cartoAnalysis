function alignFigure(currentFigure, newAlign)
% ALIGNFIGURE re-align the figure topleft or topright.
% newAlign must be 'right' or 'left'

% currentFigure = figure()
% newAlign = 'left'
% 
% alignFigure(currentFigure, newAlign)


% Get current monitor features
%monitorPositions = get(0,'MonitorPositions');
screenSize       = get(0,'ScreenSize');
screenWidth      = screenSize(3);
screenHeight     = screenSize(4);


% Let's assume any windows is windows 7
if ispc
  taskBarOffset = 40;
else
  taskBarOffset = 0;
end


switch lower(newAlign)
    case 'left'

        fromLeft   = 1;
        fromBottom = 1 + taskBarOffset;
        newWidth   = screenWidth/2;
        newHeight  = screenHeight - fromBottom;

        newPosition = [fromLeft fromBottom newWidth newHeight];

    case 'right'

        fromLeft   = 1 + screenWidth/2;
        fromBottom = 1 + taskBarOffset;
        newWidth   = screenWidth/2;
        newHeight  = screenHeight - fromBottom;

        newPosition = [fromLeft fromBottom newWidth newHeight];
    otherwise
            newPosition = get(currentFigure, 'OuterPosition');
end

set(currentFigure, 'units', 'pixels', 'OuterPosition', newPosition );


