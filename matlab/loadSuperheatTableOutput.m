function C = loadSuperheatTableOutput(filebase)

  filename = [filebase,'_ts.csv'];
    
  try
      % for output generated with a PetscViewerBinaryOpen
      A = load(filename);
      B = loadSuperheatOutput([filebase,'_0000']);
  catch
      errstr = sprintf('\nloadSuperheatTableOutput: cannot load file %s',filename);
      error(errstr); 
  end

  C.n   = A(:,1);
  C.t   = A(:,2);
  C.Cs0 = A(:,3);
  C.Cs1 = A(:,4);
  C.R   = exp(A(:,5));
  C.Cl  = A(:,6);
  C.Vl  = A(:,7);
  C.phi = C.Vl./(C.Vl + C.R.^3);
  C.F   = 1 - C.R.^3;
  C = comparisonSolutions(B,C);

  function C = comparisonSolutions(A,C)
      Stk = A.par.St/(1/A.par.K-1);
      C.Csf = Stk*lambertw(0,exp((1 - A.par.decmpr*C.t)/Stk)/Stk);
      C.Csb = 0.5*(1 - Stk - A.par.decmpr*C.t + sqrt(4*Stk + (1 - Stk - A.par.decmpr*C.t).^2));
      C.Ff  = min((-1 + C.Csf + A.par.decmpr*C.t)/A.par.St,1);
      C.Fb  = min((-1 + C.Csb + A.par.decmpr*C.t)/A.par.St,1);
  end
  
end
