function [A B C] = runSuperheatArray(fname)

  % reference parameter set
  par.decmpr = 1;
  par.K = 1e-2;
  par.St = 3;
  par.epsphi0 = 1e-3;
  par.dt = 1e-4;
  par.Fmax = 0.2;
  par.ni = 500;
  par.nout = 1e9;
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
  
  % decompression rate series
  %figure(1);
  dcr = logspace(-2,2,5);
  dt  = par.dt./dcr;
  Fm  = min(par.Fmax*(log10(dcr)+3),0.95);
  for i=1:length(dcr)
      fnme = [par.namebase,'_decmpr_',num2str(i,'%3.3d')]
      if ~exist([fnme,'_ts.csv'],'file')
          rcall = [call,' -decmpr ',num2str(dcr(i)),' -dt ',num2str(dt(i)),' -Fmax ',num2str(Fm(i)),' -filename ',fnme,' > ',fnme,'.out'];
          unix(rcall);
      end
      A{i} = loadSuperheatTableOutput(fnme);
      A_legent{i} = ['$\log\dot{\mathcal{P}}=',num2str(log10(dcr(i))),'$'];
% $$$       subplot(2,1,1); loglog(A{i}.F,A{i}.Cs0-A{i}.Cs1,'-','linewidth',2); hold on;
% $$$       subplot(2,1,2); q(i) = semilogx(A{i}.t,A{i}.F,'-','linewidth',2); hold on;
% $$$       subplot(2,1,2); semilogx(A{i}.t,A{i}.Ff,'-','linewidth',1); hold on;
% $$$       subplot(2,1,2); semilogx(A{i}.t,A{i}.Fb,'-','linewidth',1); hold on;
% $$$       legent{i} = ['$\log\dot{\mathcal{P}}=',num2str(log10(dcr(i))),'$'];
  end
% $$$   leg = legend(q, legent{:},'location','northwest');
% $$$   set(leg,'interpreter','latex');
% $$$   subplot(2,1,1); xlabel('$F$','interpreter','latex')
% $$$   ylabel('dimensionless superheating','interpreter','latex')
% $$$   set(gca,'xlim',[1e-6 1],'ylim',[1e-3 1]);
% $$$   subplot(2,1,2); xlabel('$t$','interpreter','latex')
% $$$   hold off;
    
  % decompression rate series
  %figure(2); clear legent;
  K = logspace(-5,-1,5);  
  for i=1:length(K)
      fnme = [par.namebase,'_K_',num2str(i,'%3.3d')]
      if ~exist([fnme,'_ts.csv'],'file')
          rcall = [call,' -K ',num2str(K(i)),' -filename ',fnme,' > ',fnme,'.out'];
          unix(rcall);
      end
      B{i} = loadSuperheatTableOutput(fnme);
      B_legent{i} = ['$\log K=',num2str(log10(K(i))),'$'];
% $$$       subplot(2,1,1); p(i) = loglog(B{i}.F,B{i}.Cs0-B{i}.Cs1,'-','linewidth',2); hold on;
% $$$       subplot(2,1,2); q(i) = semilogx(B{i}.t,B{i}.F,'-','linewidth',2); hold on;
% $$$       subplot(2,1,2); semilogx(B{i}.t,B{i}.Ff,'-','linewidth',1); hold on;
% $$$       subplot(2,1,2); semilogx(B{i}.t,B{i}.Fb,'-','linewidth',1); hold on;
% $$$       legent{i} = ['$\log K=',num2str(log10(K(i))),'$'];
  end
% $$$   leg = legend(q, legent{:},'location','northwest');
% $$$   set(leg,'interpreter','latex');
% $$$   subplot(2,1,1); xlabel('$F$','interpreter','latex')
% $$$   ylabel('dimensionless superheating','interpreter','latex')
% $$$   set(gca,'xlim',[1e-6 1],'ylim',[1e-3 1]);
% $$$   subplot(2,1,2); xlabel('$t$','interpreter','latex')
% $$$   hold off;
  
  % Stefan series
  %figure(3); clear legent q;
  St = 3*logspace(-1,1,3);  
  for i=1:length(St)
      fnme = [par.namebase,'_St_',num2str(i,'%3.3d')]
      if ~exist([fnme,'_ts.csv'],'file')
          rcall = [call,' -St ',num2str(St(i)),' -filename ',fnme,' > ',fnme,'.out'];
          unix(rcall);
      end
      C_legent{i} = ['St$=',num2str(St(i)),'$'];
% $$$       C{i} = loadSuperheatTableOutput(fnme);
% $$$       subplot(2,1,1); loglog(C{i}.F,C{i}.Cs0-C{i}.Cs1,'-','linewidth',2); hold on;
% $$$       subplot(2,1,2); q(i) = semilogx(C{i}.t,C{i}.F,'-','linewidth',2); hold on;
% $$$       subplot(2,1,2); semilogx(C{i}.t,C{i}.Ff,'-','linewidth',1); hold on;
% $$$       subplot(2,1,2); semilogx(C{i}.t,C{i}.Fb,'-','linewidth',1); hold on;
% $$$       legent{i} = ['St$=',num2str(St(i)),'$'];
  end
% $$$   leg = legend(q, legent{:},'location','northwest');
% $$$   set(leg,'interpreter','latex');
% $$$   subplot(2,1,1); xlabel('$F$','interpreter','latex')
% $$$   ylabel('dimensionless superheating','interpreter','latex')
% $$$   set(gca,'xlim',[1e-6 1],'ylim',[1e-3 1]);
% $$$   subplot(2,1,2); xlabel('$t$','interpreter','latex')
% $$$   hold off;
  
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
  
  fnme = [par.namebase,'_decmpr_'];
  map = colormap;
  [r c] = size(map);
  row = round(linspace(1,r,length(A)));
  colr = map(row',:);
  
  % plot of superheating versus time
  for i=1:length(A)
     p(i) = plot(A{i}.t*dcr(i),A{i}.Cs0-A{i}.Cs1,'-','linewidth',2,'color',colr(i,:)); hold on; 
  end
  xlabel('$\dot{\mathcal{P}}t$','interpreter','latex','fontsize',18);
  ylabel('dimensionless superheating','interpreter','latex','fontsize',18);
  leg = legend(p,A_legent{:},'location','southeast');
  set(leg,'interpreter','latex','fontsize',10);
  set(gca,'xlim',[0 10]);
  
  print('-dpdf',[fnme,'_SH']);
  
  set(gca,'xscale','log','yscale','log','xlim',10.^[-4 1],'ylim',10.^[-4 0.1]);
  set(leg,'location','northwest');
  
  print('-dpdf',[fnme,'_logSH']);
  
  delete(p); delete(leg);
  for i=1:length(A)
     p(i) = plot(A{i}.t*dcr(i),A{i}.F,'-','linewidth',2,'color',colr(i,:)); hold on; 
     q(i) = plot(A{i}.t*dcr(i),A{i}.Ff,'-k','linewidth',0.5); hold on; 
     qq(i) = plot(A{i}.t*dcr(i),A{i}.Fb,'--k','linewidth',0.5); hold on; 
  end
  ylabel('$F$','interpreter','latex','fontsize',18);
  leg = legend(p,A_legent{:},'location','southeast');
  set(leg,'interpreter','latex','fontsize',10);
  set(gca,'xscale','linear','yscale','linear','xlim',[0 10],'ylim',[0 1]);
  print('-dpdf',[fnme,'_F']);
  
  set(gca,'xscale','log','yscale','log','xlim',10.^[-2 1],'ylim',10.^[-4 0.1]);
  set(leg,'location','southeast');
  print('-dpdf',[fnme,'_logF']);
  
  delete(p); delete(q); delete(qq); delete(leg);
  for i=1:length(A)
     p(i) = plot(A{i}.F,A{i}.Cs0-A{i}.Cs1,'-','linewidth',2,'color',colr(i,:)); hold on; 
  end
  xlabel('$F$','interpreter','latex','fontsize',18);
  ylabel('dimensionless superheating','interpreter','latex','fontsize',18);
  leg = legend(p,A_legent{:},'location','northwest');
  set(leg,'interpreter','latex','fontsize',10);
  set(gca,'xlim',[0 1]);
  
  print('-dpdf',[fnme,'_F_SH']);
 
  
  close(f);