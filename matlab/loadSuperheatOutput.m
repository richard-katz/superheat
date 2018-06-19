function A = loadSuperheatOutput(filename)

  try
      % for output generated with a PetscViewerBinaryOpen
      A = PetscReadBinaryMatlab(filename);
  catch
      errstr = sprintf('\nloadSuperheatOutput: cannot load file %s',filename);
      error(errstr); 
  end

  A.soln.R  = exp(A.X(A.par.dofs-2));
  A.soln.Cl = A.X(A.par.dofs-1);
  A.soln.Vl = A.X(A.par.dofs);
  dr = A.soln.R/(A.par.ni-2);
  A.soln.r  = dr*([0:1:A.par.ni-1]-0.5);
  r_interp = linspace(0,A.soln.R,A.par.ni);
  A.soln.Cs = interp1(A.soln.r, A.X(1:A.par.ni), r_interp);
  A.soln.r  = r_interp;