function addSubplotIndex(f,n, h1)

a = getappdata(f, 'subplotindices');
a = unique([a,n]);
setappdata(f,'subplotindices', a);
fprintf('Adding subplot %d\n', n);
set(h1, 'LineWidth',2);