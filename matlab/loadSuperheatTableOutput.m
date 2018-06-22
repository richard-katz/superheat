function C = loadSuperheatTableOutput(filebase)

  filename = [filebase,'_ts.csv'];
    
  try
      % for output generated with a PetscViewerBinaryOpen
      B = load(filename);
  catch
      errstr = sprintf('\nloadSuperheatTableOutput: cannot load file %s',filename);
      error(errstr); 
  end

  C.n = B(:,1);
  C.t = B(:,2);
  C.Cs0 = B(:,3);
  C.Cs1 = B(:,4);
  C.lnR = B(:,5);
  C.Cl = B(:,6);
  C.Vl = B(:,7);