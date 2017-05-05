function addPointIndex(f, currentPointIdx)

a = getappdata(f, 'pointIndices');
a = unique([a,currentPointIdx]);
setappdata(f, 'pointIndices', a);
fprintf('Adding pointIdx %d\n', currentPointIdx);