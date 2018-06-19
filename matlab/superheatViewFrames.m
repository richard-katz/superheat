function A = superheatViewFrames(filebase,frames)

  for i=1:length(frames)
      filename = [filebase,sprintf('_%4.4d',frames(i))];
      A = loadSuperheatOutput(filename);
      
      p(1) = plot(A.soln.r,A.soln.Cs,'-k','linewidth',2); hold on;
  end
  hold off; axis tight; grid on;
  
  xlabel('Normalized radius, $r$','interpreter','latex');
  ylabel('Normalized concentration','interpreter','latex')
  %leg = legend(p,'$C^s$','$C^\ell$');
  %set(leg,'interpreter','latex');