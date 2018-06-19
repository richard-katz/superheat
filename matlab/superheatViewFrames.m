function A = superheatViewFrames(filebase,frames)

  yl = -eps;
  
  for i=1:length(frames)
      filename = [filebase,sprintf('_%4.4d',frames(i))];
      A = loadSuperheatOutput(filename);
      
      p(1) = plot(A.soln.r,A.soln.Cs,'-k','linewidth',2); hold on;
      yl = min(min(get(gca,'ylim')),yl);
  end
  set(gca,'ylim',[yl 0]); hold off;
  
  xlabel('Normalized radius, $\varrho$','interpreter','latex');
  ylabel('Normalized concentration','interpreter','latex')
  leg = legend(p,'$C^s$','$C^\ell$');
  set(leg,'interpreter','latex');