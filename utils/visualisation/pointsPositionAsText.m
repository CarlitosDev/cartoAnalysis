function textPosition = pointsPositionAsText(pointsPosition)


% Let's assume that the atrium is quite round. Get the centroid and
% the direction of a vector from the centroid to every point. Then,
% shift the text by 'delta' pointing outwards.
deltaText      = 2.5;

pointsCentroid = mean(pointsPosition, 1);
pointsCentered = ...
  bsxfun(@minus, pointsPosition, pointsCentroid);

vectorDirection = sign(pointsCentered);

textPosition = pointsPosition + deltaText.*vectorDirection;