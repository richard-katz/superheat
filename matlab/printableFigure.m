function f = printableFigure(width, height, units)

  if nargin==2; units = 'inches'; end
    f = figure;
    set(f,'Units',units,'Position',[0.7 12 fw fh]);
    set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(f,'Color','w','InvertHardcopy','off', 'Resize','off');
  