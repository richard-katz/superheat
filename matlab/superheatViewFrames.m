function A = superheatViewFrames(filebase,frames)

  FS = {'FontSize',20};
    
  for i=1:length(frames)
      filename = [filebase,sprintf('_%4.4d',frames(i))];
      try A(i) = loadSuperheatOutput(filename);
      catch break; end
  end
  
  C = loadSuperheatTableOutput(filebase);
  
  subplot(2,1,1);
  p(1) = plot(C.t,C.Cs0-C.Cs1,'-k','linewidth',2); hold on;
  p(2) = plot(C.t,exp(C.lnR),'-r','linewidth',2); 
  p(3) = plot(C.t,C.Cl,'-b','linewidth',2); 
  p(4) = plot(C.t,C.Vl./(C.Vl + exp(C.lnR).^3),'-m','linewidth',2);
  p(5) = plot(C.t,1-exp(C.lnR).^3,'-g','linewidth',2);
  hold off;
  ylabel('$\Delta T, R, C^\ell, \phi, F$','interpreter','latex',FS{:})
  xlabel('Dimensionless time, $t$','interpreter','latex',FS{:});
  set(gca,'ylim',[0 1]);
  
  tmin = A(1).par.t;
  tmax = A(end).par.t;
  nmax = A(end).par.n;
  tspan = tmax-tmin;
  for i=1:length(A); 
      t(i) = A(i).par.t/tspan; 
      tind(i) = find(C.n==A(i).par.n);
  end
  
  map = colormap;
  [r c] = size(map);
  row = max(round(t*r),1);
  colr = map(row',:);
  
  for i=1:length(A)
      subplot(2,1,2);
      plot(A(i).soln.r,A(i).soln.Cs,'-k','linewidth',2,'color',colr(i,:)); hold on;
      %subplot(3,1,3);
      %q(i) = plot(A(i).soln.r,A(i).soln.Cs-A(i).soln.Cs(end),'-k','linewidth',2,'color',colr(i,:)); hold on;
  end
  hold off; axis tight; grid on; 
  
  subplot(2,1,1);  hold on;
  scatter(C.t(tind),C.Cs0(tind)-C.Cs1(tind),[80],colr,'linewidth',2); hold off
  leg = legend(p,'$\Delta T(t)$','$R(t)$','$C^\ell(t)$','$\phi(t)$','$F$');
  set(leg,'interpreter','latex',FS{:});
  
  ti = ['$\dot{\mathcal{P}}=$',num2str(-A(1).par.decmpr),'$,\;K=$',...
        num2str(A(1).par.K,'%.1e'),', St$=$',num2str(A(1).par.St),', $\phi\le$',num2str(A(1).par.phi)];
  title(ti,'interpreter','latex',FS{:});
  
  subplot(2,1,2)
  ylabel('Normalized concentration','interpreter','latex',FS{:})
  xlabel('Normalized radius, $r$','interpreter','latex',FS{:});

  %figure;
  %semilogy(C.t,C.Cl,'-b','linewidth',2);
  
  %ylabel('Dimensionless superheating','interpreter','latex')
  %xlabel('Normalized radius, $r$','interpreter','latex');
  %leg = legend(p,'$C^s$','$C^\ell$');
  %set(leg,'interpreter','latex');