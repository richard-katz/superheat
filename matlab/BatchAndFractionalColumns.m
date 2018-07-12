function A = BatchAndFractionalColumns(par)
    
    if nargin==0
        par.D = 0.01;           % distribution coefficient
        par.c = 1200;           % heat capacity J/kg/K
        par.alpha = 3e-5;       % expansivity /K
        par.rho = 3000;         % density
        par.g = 10;             % gravity
        par.L = 5e5;            % latent heat J/kg
        par.C = 6.5e6;          % clapeyron Pa/K
        par.M = -4;             % solidus slope, K per ppm volatile
        par.c0 = 100;           % volatile concentration, ppm
        par.TsP0 = 1100 + 273;  % solidus at P=0
        par.Ts0 = 1300+273;     % mantle temperature
        par.phi0 = 0.01;        % reference porosity
        par.epsilon = 0.01;     % melt/solid speed ratio
        par.n = 2;              % permeability exponent
    end
    
    Mc0 = par.M*par.c0;
    par.z0 = (par.Ts0 - par.TsP0 - Mc0)/(par.rho*par.g/par.C);
    
    A.z = linspace(0,par.z0,300);
    Tcc = par.Ts0 - par.rho*par.g/par.C*A.z;
    
    % batch
    A1 = par.L*par.c0/par.c/par.M/(1-1/par.D);
    A2 = -par.g/par.M*(par.alpha*par.Ts0/par.c - par.rho/par.C);
    A.cs = 0.5*(par.c0 - A1/par.c0 + A2*A.z + sqrt(4*A1 + (A2*A.z + par.c0 - A1/par.c0).^2));
    A.F = par.c/par.L*(par.rho*par.g*A.z/par.C + ...
                       par.M*(par.c0-A.cs) - par.alpha*par.g*par.Ts0*A.z/par.c);
    A.T = par.Ts0 - par.rho*par.g*A.z/par.C + par.M*(A.cs - par.c0);
    A.phi = par.epsilon/2*(sqrt(1 + 4*A.F/par.epsilon/par.phi0) - 1);
    
    % fractional
    A1 = A1/par.c0;
    A.csf = A1*lambertw(0,(par.c0*exp((par.c0 + A2*A.z)/A1))/A1);
    A.Ff = par.c/par.L*(par.rho*par.g*A.z/par.C + ...
                        par.M*(par.c0-A.csf) - par.alpha*par.g*par.Ts0*A.z/par.c);
    A.Tf = par.Ts0 - par.rho*par.g*A.z/par.C + par.M*(A.csf-par.c0);
    A.phif = par.epsilon/2*(sqrt(1 + 4*A.Ff/par.epsilon/par.phi0) - 1);
  
    A.par = par;
    
end
