function saveFigureAsPNG(figTitle, figPath)


f = gcf();
set(f,'units','normalized','outerposition',[0 0 1 1]);

% save figure
set(f, 'Name' , figTitle);


figPathName = fullfile(figPath, [figTitle, '.png']);
hgexport(f, figPathName, hgexport('factorystyle'), 'Format', 'png');
close(f);


