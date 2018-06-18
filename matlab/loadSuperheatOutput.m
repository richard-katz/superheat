function A = loadSuperheatOutput(filename)

  try
      % for output generated with a PetscViewerBinaryOpen
      A = PetscReadBinaryMatlab(filename);
  catch
      errstr = sprintf('\nloadSuperheatOutput: cannot load file %s',filename);
      error(errstr); 
  end

  A.soln.Cs = A.X(1:A.par.ni);
  A.soln.R  = exp(A.X(A.par.ni+1));
  dr = A.soln.R/(A.par.ni-2);
  A.soln.r  = dr*([0:1:A.par.ni-1]-0.5);
