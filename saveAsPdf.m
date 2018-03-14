function saveAsPdf(fig_number, filename)

% saves the figure(fig_number) as a pdf (filename.pdf)
% 
% fig_number    The number of the figure to save.
% 
% filename      String with the name for the file, extension can be ommited
%               it will be added automatically.
% 
% Author: Rodrigo Abrajan

h = figure(fig_number);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if(strfind(filename,'.pdf')==length(filename)-3)
    print(h,filename,'-dpdf','-r0')
else
    print(h,strcat(filename,'.pdf'),'-dpdf','-r0')
end