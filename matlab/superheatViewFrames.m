function A = superheatViewFrames(filebase,frames)

  for i=1:length(frames)
      filename = [filebase,sprintf('_%4.4d',frames(i))];
      A(i) = loadSuperheatOutput(filename);
  end
  
  tmin = A(1).par.t;
  tmax = A(end).par.t;
  tspan = tmax-tmin;
  for i=1:length(A); t(i) = A(i).par.t/tspan; end
  
  map = colormap;
  [r c] = size(map);
  row = max(round(t*r),1);
  colr = map(row',:);
  
  for i=1:length(A)
      subplot(2,1,1);
      p(i) = plot(A(i).soln.r,A(i).soln.Cs,'-k','linewidth',2,'color',colr(i,:)); hold on;
      subplot(2,1,2);
      q(i) = plot(A(i).soln.r,A(i).soln.Cs-A(i).soln.Cs(end),'-k','linewidth',2,'color',colr(i,:)); hold on;
  end
  hold off; axis tight; grid on; 
  
  subplot(2,1,1)
  ylabel('Normalized concentration','interpreter','latex')
  subplot(2,1,2)
  ylabel('Dimensionless superheating','interpreter','latex')
  xlabel('Normalized radius, $r$','interpreter','latex');
  %leg = legend(p,'$C^s$','$C^\ell$');
  %set(leg,'interpreter','latex');