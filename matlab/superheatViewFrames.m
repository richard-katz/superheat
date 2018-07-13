function [A B C] = superheatViewFrames(filebase,frames)

  FS = {'FontSize',20};
    
  for i=1:length(frames)
      filename = [filebase,sprintf('_%4.4d',frames(i))];
      try A(i) = loadSuperheatOutput(filename);
      catch display(['failed to load file ',filename]);  break; 
      end
  end
  
  C = loadSuperheatTableOutput(filebase);
  B = comparisonSolutions(A(1),C);
  
  %subplot(2,1,1); plot(C.t,B.Csf,C.t,B.Csb);
  %subplot(2,1,2); plot(C.t,B.Ff,C.t,B.Fb); return

  
  subplot(3,1,1);
  p(1) = plot(C.t,C.Cs0-C.Cs1,'-k','linewidth',2); hold on;
  p(2) = plot(C.t,exp(C.lnR),'-r','linewidth',2); 
  p(3) = plot(C.t,C.Cl,'-b','linewidth',2); 
  p(4) = plot(C.t,C.Vl./(C.Vl + exp(C.lnR).^3),'-m','linewidth',2);
  hold off;
  ylabel('$\Delta T, R, C^\ell, \phi$','interpreter','latex',FS{:})
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
      subplot(3,1,3);
      plot(A(i).soln.r,A(i).soln.Cs,'-k','linewidth',2,'color',colr(i,:)); hold on;
  end
  hold off; axis tight; grid on; 
  
  subplot(3,1,1);  hold on;
  scatter(C.t(tind),C.Cs0(tind)-C.Cs1(tind),[80],colr,'linewidth',2); hold off
  leg = legend(p,'$\Delta T(t)$','$R(t)$','$C^\ell(t)$','$\phi(t)$');
  set(leg,'interpreter','latex',FS{:});
  
  ti = ['$\dot{\mathcal{P}}=$',num2str(-A(1).par.decmpr),'$,\;K=$',...
        num2str(A(1).par.K,'%.1e'),', St$=$',num2str(A(1).par.St),', $\epsilon\phi_0=$',num2str(A(1).par.epsphi0)];
  title(ti,'interpreter','latex',FS{:});

  subplot(3,1,2); hold on;
  p(1) = plot(C.t,1-exp(C.lnR).^3,'-b','linewidth',2);
  p(3) = plot(C.t,B.Fb,'-k','linewidth',2);
  p(2) = plot(C.t,B.Ff,'--r','linewidth',2);
  p(4) = plot(C.t,C.Vl./(C.Vl + exp(C.lnR).^3),'-m','linewidth',2);
  leg = legend(p(1:4),'grain $F(t)$','fractional $F(t)$','batch $F(t)$','$\phi(t)$');
  set(leg,'interpreter','latex',FS{:},'location','northwest');
  
  
  subplot(3,1,3)
  ylabel('Normalized concentration','interpreter','latex',FS{:})
  xlabel('Normalized radius, $r$','interpreter','latex',FS{:});

  function B = comparisonSolutions(A,C)
      Stk = A.par.St/(1/A.par.K-1);
      B.Csf = Stk*lambertw(0,exp((1 - A.par.decmpr*C.t)/Stk)/Stk);
      B.Csb = 0.5*(1 - Stk - A.par.decmpr*C.t + sqrt(4*Stk + (1 - Stk - A.par.decmpr*C.t).^2))
      B.Ff  = min((-1 + B.Csf + A.par.decmpr*C.t)/A.par.St,1);
      B.Fb  = min((-1 + B.Csb + A.par.decmpr*C.t)/A.par.St,1);
  end
  
end