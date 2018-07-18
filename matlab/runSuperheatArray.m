function A = runSuperheatArray(fname)

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
  figure(1);
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
      loglog(A{i}.F,A{i}.Cs0-A{i}.Cs1,'-','linewidth',2); hold on;
  end
  set(gca,'xlim',[1e-6 1],'ylim',[1e-3 1]);
  hold off;
    
  % decompression rate series
  figure(2);
  K = logspace(-5,-1,5);  
  for i=1:length(K)
      fnme = [par.namebase,'_K_',num2str(i,'%3.3d')]
      if ~exist([fnme,'.csv'],'file')
          rcall = [call,' -K ',num2str(K(i)),' -filename ',fnme,' > ',fnme,'.out'];
          unix(rcall);
      end
      B{i} = loadSuperheatTableOutput(fnme);
      loglog(A{i}.F,A{i}.Cs0-A{i}.Cs1,'-','linewidth',2); hold on;
  end
  set(gca,'xlim',[1e-6 1],'ylim',[1e-3 1]);
  hold off;