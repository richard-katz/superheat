function [A B C] = runSuperheat2DArray(fname)

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
  opts = [' -St ',num2str(par.St)];
  opts = [opts,' -epsphi0 ',num2str(par.epsphi0)];
  opts = [opts,' -Fmax ',num2str(par.Fmax)];
  opts = [opts,' -dt ',num2str(par.dt)];
  opts = [opts,' -ni ',num2str(par.ni)];
  opts = [opts,' -nout ',num2str(par.nout)];
  call = [prog,opts];
  
  % decompression rate series
  %figure(1);
  dcr = logspace(-2,2,2);
  K   = logspace(-4,-0.5,2);
  dt  = par.dt./dcr;
  Fm  = min(par.Fmax*(log10(dcr)+3),0.95)
  if ~exist('sh2darray.mat','file')
      for i=1:length(dcr)
          for j=1:length(K)
              fnme = [par.namebase,'_decmpr_',num2str(i,'%3.3d'),'_K_',num2str(j,'%3.3d')]
              if ~exist([fnme,'_ts.csv'],'file')
                  rcall = [call,' -decmpr ',num2str(dcr(i)),' -dt ',num2str(dt(i)) ... 
                           ,' -K ',num2str(K(j)) ...
                           ,' -Fmax ',num2str(Fm(i)),' -filename ',fnme,' > ',fnme,'.out'];
                  unix(rcall);
              end
              A{i,j} = loadSuperheatTableOutput(fnme);
          end
      end    
      save sh2darray A;
  end
  
  load sh2darray;
  for i=1:length(dcr)
      for j=1:length(K)
          plot(A{i,j}.Fb,A{i,j}.Cs0-A{i,j}.Cs1); hold on;
          mash(i,j) = max(A{i,j}.Cs0 - A{i,j}.Cs1)
      end
  end
  %imagesc(log(dcr),log(K),mash);
    

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
  
