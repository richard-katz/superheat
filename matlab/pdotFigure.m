clear; close all;

  % axis dimensions, inches
  axh = 4;    % axis height
  axw = 4;    % axis width
  axb = 0.7;  % axis bottom spacing
  axt = 0.1;  % axis top spacing
  axl = 0.9;  % axis left spacing
  axr = 0.8;  % axis right spacing
  fw = axl + axw + axr;
  fh = axb + axh + axt;
  f = printableFigure(fw,fh);
  ax = axes('units','inches','position',[axl axb axw axh]);
  
  par.rho = 3300;   % kg/m^3
  par.g = 10;       % m/sec^2
  par.D = 1e-18;    % m^2/sec
  par.Mc0s = 200;   % K
  par.clap = 6.5e6; % Pa/K
  
  W0 = logspace(0,3,100); % cm/year
  R0 = logspace(-1,2,100); % mm
  
  W0 = W0/100/pi/1e7; % m/sec
  R0 = R0/1e3;        % m
  
  [W R] = meshgrid(W0,R0);
  
  Pdot = par.rho*par.g*W.*R.^2/(par.D*par.clap*par.Mc0s);
  pdmin = min(min(Pdot));
  pdmax = max(max(Pdot));
  clevs = logspace(log10(pdmin),log10(pdmax),12);
  
  contourf(W0*pi*1e9,R0*1e3,log10(Pdot),[-5:5]); hold on;
  contour(W0*pi*1e9,R0*1e3,log10(Pdot),[0 0],'color','k','linewidth',2);
  set(gca,'yscale','log','xscale','log','xtick',10.^[-1:3],...
          'xticklabel',10.^[-1:3],'ytick',10.^[-2:2],'yticklabel',10.^[-2:2]);
  xlabel('$W_0$, cm/yr','interpreter','latex','fontsize',18)
  ylabel('$R_0$, mm','interpreter','latex','fontsize',18)

  
  cb = colorbar('location','East','units','inches');
  cpos = get(cb,'position');
  cpos(1) = axl + axw + 0.3;
  set(cb,'position',cpos);
  
  print -dpdf pdotFigure
  close all;
  !open pdotFigure.pdf