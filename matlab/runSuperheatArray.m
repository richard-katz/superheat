function A = runSuperheatArray(N, fname)

  par.N = N;
  par.phi = 0;
  par.decmpr = logspace(-1,2,par.N);
  par.dt = 1e-4;
  par.K = 1e-2;
  par.St = 10;
  par.nout = 1e9;
  par.namebase = fname;
  
  for i=1:par.N
      fnme = [par.namebase,'_',num2str(i,'%3.3d')]
      call = ['./superheat -snes_fd -dt ',num2str(par.dt)];
      call = [call,' -phi ',num2str(par.phi)];
      call = [call,' -decmpr ',num2str(par.decmpr(i))];
      call = [call,' -K ',num2str(par.K)];
      call = [call,' -St ',num2str(par.St)];
      call = [call,' -nout ',num2str(par.nout)];
      call = [call,' -filename ',fnme,' > ',fnme,'.out'];
      unix(call);
      A{i} = loadSuperheatTableOutput(fnme);
  end