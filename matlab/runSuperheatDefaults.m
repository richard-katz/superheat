function A = runSuperheatDefaults(fname)

  % reference parameter set
  par.decmpr = 10;
  par.K = 1e-2;
  par.St = 3;
  par.epsphi0 = 1e-3;
  par.dt = 1e-4;
  par.Fmax = 0.50;
  par.ni = 500;
  par.nout = 200;
  par.namebase = fname;
  prog = '../src/superheat -snes_fd';
  opts = [' -decmpr ',num2str(par.decmpr)];
  opts = [opts,' -K ',num2str(par.K)];
  opts = [opts,' -St ',num2str(par.St)];
  opts = [opts,' -epsphi0 ',num2str(par.epsphi0)];
  opts = [opts,' -Fmax ',num2str(par.Fmax)];
  opts = [opts,' -dt ',num2str(par.dt)];
  opts = [opts,' -ni ',num2str(par.ni)];
  opts = [opts,' -nout ',num2str(par.nout)];
  call = [prog,opts];

  
  fnme = [par.namebase,'_defaults'];
  if ~exist([fnme,'_ts.csv'],'file')
      unix([call,' -filename ',fnme,' > ',fnme,'.out']);
  end
  
  C = loadSuperheatTableOutput(fnme);
  
  for j=1:1000
      filename = [fnme,sprintf('_%4.4d',j)];
      try A(j) = loadSuperheatOutput(filename);
      catch break; end
  end  
  
  tmin = A(1).par.t;
  tmax = A(end).par.t;
  nmax = A(end).par.n;
  tspan = tmax-tmin;
  for i=1:length(A); 
      t(i) = (A(i).par.t-tmin)/tspan; 
      tind(i) = find(C.n==A(i).par.n);
  end
  
  
  %%%%%%%%%% MAKE PLOTS %%%%%%%%%%%%%
  
  % axis dimensions, inches
  axh = 3;    % axis height
  axw = 5;    % axis width
  axb = 0.6;  % axis bottom spacing
  axt = 0.1;  % axis top spacing
  axl = 0.8;  % axis left spacing
  axr = 0.2;  % axis right spacing
  fw = axl + axw + axr;
  fh = axb + axh + axt;
  f = printableFigure(fw,fh);
  
  ax = axes('units','inches','position',[axl axb axw axh]);
  
  map = colormap;
  [r c] = size(map);
  row = max(round(t*r),1);
  colr = map(row',:);
  
  for i=1:length(A)
      plot(A(i).soln.r,A(i).soln.Cs,'-k','linewidth',2,'color',colr(i,:)); hold on;
  end
  hold off; axis tight; grid on; 
  ylabel('$C^s$','interpreter','latex','fontsize',18)
  xlabel('$r$','interpreter','latex','fontsize',18)
  text(0.02, 0.02,'$\dot{\mathcal{P}}=10$','interpreter','latex','fontsize',18,...
       'horizontalalignment','left','verticalalignment','bottom','backgroundcolor','w','units','normalized')
  
  print('-dpdf',[fnme,'_r_Cst']);
  
  set(gca,'yscale','log');
  
  print('-dpdf',[fnme,'_r_logCst']);
  
  clf;  ax = axes('units','inches','position',[axl axb axw 0.9*axh]);

  p(1) = plot(C.t,C.Cs0-C.Cs1,'-k','linewidth',2); hold on;
  p(3) = plot(C.t,C.Cl,'-b','linewidth',2); 
  p(2) = plot(C.t,C.R,'-r','linewidth',2); 
  p(4) = plot(C.t,C.F,'-m','linewidth',2);
  scatter(C.t(tind),C.Cs0(tind)-C.Cs1(tind),[80],colr,'linewidth',2); hold off
  leg = legend(p,'$\Delta T(t)$','$C^\ell(t)$','$R(t)$','$F(t)$','orientation','horizontal');
  set(leg,'interpreter','latex','fontsize',14,'units','inches');
  lpos = get(leg,'position'); lpos(2) = axb + 0.9*axh + 0.05; set(leg,'position',lpos);
  hold off;
  ylabel('$\Delta T, C^\ell,  R, F$','interpreter','latex','fontsize',18)
  xlabel('Dimensionless time, $t$','interpreter','latex','fontsize',18);
  set(gca,'ylim',[0 1]);
  
  
  print('-dpdf',[fnme,'_t_various']);
  
  close all;