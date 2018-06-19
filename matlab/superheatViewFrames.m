function A = superheatViewFrames(filebase,frames)

  filename = [filebase,sprintf('_%4.4d',frames(1))];
  A = loadSuperheatOutput(filename);
  
  p(1) = plot(A.soln.r,A.soln.Cs,'-k','linewidth',2);
  yl = get(gca,'ylim');
  set(gca,'ylim',[yl(1) 0]);
  
  xlabel('Normalized radius, $\varrho$','interpreter','latex');
  ylabel('Normalized concentration','interpreter','latex')
  leg = legend(p,'$C^s$','$C^\ell$');
  set(leg,'interpreter','latex');