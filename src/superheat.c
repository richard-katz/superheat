#include "superheat.h"

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
/*-----------------------------------------------------------------------*/
{
  AppCtx         user;  
  PetscErrorCode ierr;
  
  /* initialize */
  PetscInitialize(&argc,&argv,(char*)0,PETSC_NULL);
  ierr = SetUpParameters(&user);CHKERRQ(ierr);
  ierr = PetscOptionsInsert(NULL,&argc,&argv,PETSC_NULL);CHKERRQ(ierr);
  ierr = SetUpDataStructures(&user);CHKERRQ(ierr);

  /* form IC, solve and output */
  ierr = SetUpInitialGuess(&user);CHKERRQ(ierr);
  ierr = WriteOutput(&user);CHKERRQ(ierr);
  ierr = DoSolve(&user);CHKERRQ(ierr);

  /* clean up */
  ierr = CleanUpDataStructures(&user);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormResidual"
PetscErrorCode FormResidual(SNES snes, Vec X, Vec Res, void *ptr)
/* ------------------------------------------------------------------- */
{
  PetscErrorCode  ierr;
  AppCtx          *user = (AppCtx*)ptr;
  Parameter       *par = user->param;
  PetscReal const *x, *xo;
  PetscReal       dt = user->param->dt;
  PetscReal       dr, r, R, *res, Rdot, Cdot, Vdot, Rsq, Rcu;
  PetscReal       CsdotR, Cl, Cldot, Vl, CsR, GradCsR, phi;
  PetscInt        i, iR, iCl, iV, iCs;
  
  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(X,&x); CHKERRQ(ierr);
  ierr = VecGetArrayRead(user->Xo,&xo); CHKERRQ(ierr);
  ierr = VecGetArray(Res,&res); CHKERRQ(ierr);
  dr = 1./(user->param->ni-2);

  /* extract ODE variables */
  iR  = user->param->dofs - N_ODES;     R   = exp(x[iR]);  
  iCl = user->param->dofs - N_ODES + 1; Cl  = x[iCl]; 
  iV  = user->param->dofs - N_ODES + 2; Vl  = x[iV]; 
  iCs = user->param->ni   - 1;          CsR = 0.5*(x[iCs] + x[iCs-1]);

  /* other useful quantities */
  GradCsR = (x[iCs] - x[iCs-1])/dr;
  Rdot    = (x[iR] - xo[iR])/dt;
  CsdotR  = (0.5*(x[iCs]+x[iCs-1]) - 0.5*(xo[iCs]+xo[iCs-1]))/dt;
  Vdot    = (x[iV] - xo[iV])/dt;
  Cldot   = (x[iCl] - xo[iCl])/dt;
  Rsq = R*R; Rcu = Rsq*R;
  
  /* radius ODE (complete) */
  res[iR] = Rdot - (-par->decmpr - CsdotR)
          / (3*par->St/(1 + Vl/pow(R,3)) - GradCsR);
  
  /* liquid volume ODE (complete) */
  phi = 0.5*par->epsphi0*(sqrt(1 + 4*(1-Rcu)/par->epsphi0) - 1);
  res[iV] = Vl - Rcu*phi/(1-phi);
  
  /* liquid concentration ODE (complete) */
  res[iCl] = Vl*Cldot/Rcu/3
           + Rdot*(par->K-1)*Cl
           + par->K*GradCsR/Rsq;

  /* diffusion boundary conditions (complete) */
  res[0]    = x[0] - x[1];
  res[iCs]  = CsR - Cl;

  /* diffusion PDE (complete) */
  for (i=1; i<iCs; i++) {
    r = dr*(i-0.5);  Cdot = (x[i] - xo[i])/dt;
    res[i] = -Cdot + r*Rdot*(x[i+1] - x[i-1])/(2*dr) +
           + ((x[i-1] - 2*x[i] + x[i+1])/dr/dr + (x[i+1] - x[i-1])/r/dr)/Rsq;
  }
  
  ierr = VecRestoreArray(Res,&res); CHKERRQ(ierr);
  ierr = VecScale(Res,dr);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(user->Xo,&xo); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(X,&x); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat B, void *ptr)
/* ------------------------------------------------------------------- */
{
  PetscErrorCode  ierr;
  AppCtx          *user = (AppCtx*)ptr;
  Parameter       *par = user->param;
  PetscReal const *x, *xo;
  PetscInt        row, col[5], is, ie, iR, iCl, iV, iCs;
  PetscReal       A[5], dr, r, dt, R, Rdot, delsq, Cl, Cldot, Vdot, Vl, CdotR, Cs;
  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(X,&x); CHKERRQ(ierr);
  ierr = VecGetArrayRead(user->Xo,&xo); CHKERRQ(ierr);

  is = 0; ie = user->param->ni-1;
  dr = 1./(user->param->ni-2);
  dt = user->param->dt;

  /* extract ODE variables */
  iCs = user->param->ni - 1;            Cs = 0.5*(x[iCs] + x[iCs-1]);
  iR  = user->param->dofs - N_ODES;     R  = exp(x[iR]);  
  iCl = user->param->dofs - N_ODES + 1; Cl = x[iCl]; 
  iV  = user->param->dofs - N_ODES + 2; Vl = x[iV]; 

  /* radius ODE (complete) */
  R = exp(x[iR]);  Rdot = (x[iR] - xo[iR])/dt; 
  CdotR = (Cs - 0.5*(xo[iCs]+xo[iCs-1]))/dt;
  /*R */       col[0] = iR;    A[0] = 1/dt - 3*pow(R,3)*par->St*Vl*(par->decmpr + CdotR)
			  	    / pow(pow(R,3)*((x[iCs] - x[iCs-1])/dr - par->St)
					  + (x[iCs] - x[iCs-1])/dr*Vl, 2);
  /*Cs left*/  col[1] = ie-1;  A[1] = +dr*(xo[iCs] + xo[iCs-1] - 2*(x[iCs]   + par->decmpr*dt) + par->St*dr)
				    / (2*dt*pow(par->St*dr - x[iCs] + x[iCs-1],2));
  /*Cs right*/ col[2] = ie;    A[2] = -dr*(xo[iCs] + xo[iCs-1] - 2*(x[iCs-1] + par->decmpr*dt) - par->St*dr)
				    / (2*dt*pow(par->St*dr - x[iCs] + x[iCs-1],2));
  /*Vl */      col[3] = iV;    A[3] = par->St*(par->decmpr + CdotR)/pow(R,3)
				    / pow((x[iCs] - x[iCs-1])/dr*(1 + Vl/pow(R,3)) - par->St, 2);
  ierr = MatSetValues(J,1,&iR,4,col,A,INSERT_VALUES);CHKERRQ(ierr);

  /* liquid volume ODE (BROKEN) */
  Vl = x[iV];  Vdot = (x[iV] - xo[iV])/dt;
  /* Vl */      col[0] = iV;   A[0] = 1/dt;
  /* R  */      col[1] = iR;   A[1] = 3*pow(R,3)*(3*Rdot + 1/dt);
  ierr = MatSetValues(J,1,&iV,2,col,A,INSERT_VALUES);CHKERRQ(ierr);
  
  /* liquid concentration ODE (complete) */
  Cl = x[iCl]; Cldot = (x[iCl] - xo[iCl])/dt;
  /*Cl */      col[0] = iCl;  A[0] = Vl/dt/pow(R,3)/3 + Rdot*(par->K-1); 
  /*R  */      col[1] = iR;   A[1] = - Vl*Cldot/pow(R,3)
				     + (par->K-1)*Cl/dt
				     - 2*par->K*(x[ie]-x[ie-1])/dr/R/R;
  /*Vl */      col[2] = iV;   A[2] = Cldot/pow(R,3)/3;
  /*Cs left*/  col[3] = ie-1; A[3] = -par->K/dr/R/R;
  /*Cs right*/ col[4] = ie;   A[4] = +par->K/dr/R/R;
  ierr = MatSetValues(J,1,&iCl,5,col,A,INSERT_VALUES);CHKERRQ(ierr);

  /* diffusion boundary conditions (complete) */
  col[0]=is;   A[0] = +1;
  col[1]=is+1; A[1] = -1;
  ierr = MatSetValues(J,1,&is,2,col,A,INSERT_VALUES);CHKERRQ(ierr);
  col[0]=ie-1; A[0] = +0.5;
  col[1]=ie;   A[1] = +0.5;
  col[2]=iCl;  A[2] = -1;
  ierr = MatSetValues(J,1,&ie,3,col,A,INSERT_VALUES);CHKERRQ(ierr);
  
  /* diffusion PDE (complete) */
  for (row=is+1; row<ie; row++) {
    r = dr*(row-0.5);
    col[0] = row-1; col[1] = row; col[2] = row+1; col[3] = iR;
    /*left*/ A[0] = (1/dr - 1/r)/dr/R/R - r*Rdot/(2*dr);
    /*cntr*/ A[1] = -2/dr/dr/R/R - 1/dt;
    /*rght*/ A[2] = (1/dr + 1/r)/dr/R/R + r*Rdot/(2*dr);
    delsq = ((x[row-1] - 2*x[row] + x[row+1])/dr + (x[row+1] - x[row-1])/r)/dr/R/R;    
    /*R   */ A[3] = r*(x[row+1] - x[row-1])/(2*dr)/dt - 2*delsq;
    ierr = MatSetValues(J,1,&row,4,col,A,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatScale(J,dr);CHKERRQ(ierr);
  
  ierr = VecRestoreArrayRead(user->Xo,&xo); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(X,&x); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "SetUpParameters"
PetscErrorCode SetUpParameters(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  Parameter      *par;
  PetscFunctionBeginUser;

  /* Set up data structures */
  user->comm = PETSC_COMM_WORLD;
  ierr = PetscBagCreate(user->comm,sizeof(Parameter),&(user->bag)); CHKERRQ(ierr); 
  ierr = PetscBagGetData(user->bag,(void**)&user->param);CHKERRQ(ierr);
  ierr = PetscBagSetName(user->bag,"par","(parameters for crystal superheating problem)");CHKERRQ(ierr);
  par = user->param;

  /* Register numerical model parameters */
  ierr = PetscBagRegisterInt(user->bag,&par->ni,500,"ni","Number of grid points"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->dofs,par->ni+N_ODES,"dofs","<DO NOT SET> Total number of degrees of freedom"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->ns,1e5,"ns","Number of time steps"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->nout,100,"nout","Output step interval"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->n,0,"n","<DO NOT SET> Current time-step number"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->N,0,"N","<DO NOT SET> Current output frame number"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->t,0,"t","<DO NOT SET> Time");CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->tmax,1e9,"tmax","Maximum time");CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->Fmax,0.25,"Fmax","Maximum degree of melting");CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->dt,1e-3,"dt","Time-step size");CHKERRQ(ierr);
  ierr = PetscBagRegisterString(user->bag,&par->filename,FNAME_LENGTH,"test","filename","Name of output file");CHKERRQ(ierr);

  /* Register physical parameters */
  ierr = PetscBagRegisterReal(user->bag,&par->K,1e-2,"K","Parition coefficient");CHKERRQ(ierr);  
  ierr = PetscBagRegisterReal(user->bag,&par->decmpr,1,"decmpr","Dimensionless decompression rate");CHKERRQ(ierr);  
  ierr = PetscBagRegisterReal(user->bag,&par->St,3,"St","Stefan number");CHKERRQ(ierr);  
  ierr = PetscBagRegisterReal(user->bag,&par->epsphi0,1e-3,"epsphi0","Velocity ratio W0/w0 times reference porosity");CHKERRQ(ierr);  
  
  /* Display parameters */
  ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscBagView(user->bag,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  
  PetscOptionsSetValue(NULL,"-snes_monitor","");
  PetscOptionsSetValue(NULL,"-snes_converged_reason","");
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "WriteOutput"
PetscErrorCode WriteOutput(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  char           filename[FNAME_LENGTH];
  PetscViewer    viewer;
  PetscFunctionBeginUser;

  /* create frame output */
  sprintf(filename,"%s_%4.4d",user->param->filename,user->param->N++);
  ierr = PetscPrintf(user->comm,"Generating output file: \"%s\"\n",filename);
  ierr = PetscViewerBinaryOpen(user->comm,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = PetscBagView(user->bag,viewer);CHKERRQ(ierr); 
  ierr = VecView(user->X,viewer);CHKERRQ(ierr);
  ierr = VecView(user->R,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "SetUpInitialGuess"
PetscErrorCode SetUpInitialGuess(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  PetscInt       ii[N_ODES];
  PetscReal      vals[N_ODES];
  PetscFunctionBeginUser;
  ierr = VecSet(user->X,1);CHKERRQ(ierr);
  /* ln R */ ii[0] = user->param->dofs - N_ODES; vals[0] = 0;
  /* Cl   */ ii[1] = ii[0]+1;                    vals[1] = 1;
  /* Vl   */ ii[2] = ii[0]+2;                    vals[2] = 0;
  ierr = VecSetValues(user->X,N_ODES,ii,vals,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecCopy(user->X,user->Xo);CHKERRQ(ierr);
  ierr = TimestepTableAddEntry(user);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "SetUpDataStructures"
PetscErrorCode SetUpDataStructures(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  Parameter      *par = user->param;
  char           filename[FNAME_LENGTH];
  PetscFunctionBeginUser;

  /* vectors */
  ierr = VecCreate(user->comm,&user->X);CHKERRQ(ierr);
  ierr = VecSetType(user->X,VECSEQ);CHKERRQ(ierr);
  ierr = VecSetSizes(user->X,PETSC_DECIDE,par->dofs);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user->X,"X"); CHKERRQ(ierr);
  ierr = VecSetFromOptions(user->X);CHKERRQ(ierr);
  ierr = VecSetUp(user->X);CHKERRQ(ierr);
  ierr = VecDuplicate(user->X,&user->Xo);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user->Xo,"Xo"); CHKERRQ(ierr);
  ierr = VecDuplicate(user->X,&user->R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user->R,"R"); CHKERRQ(ierr);

  /* matricies */
  ierr = MatCreate(user->comm,&user->J);CHKERRQ(ierr);
  ierr = MatSetSizes(user->J,PETSC_DECIDE,PETSC_DECIDE,par->dofs,par->dofs);CHKERRQ(ierr);
  ierr = MatSetType(user->J,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(user->J,5,NULL); CHKERRQ(ierr);
  ierr = MatSetFromOptions(user->J);CHKERRQ(ierr);
  ierr = MatSetUp(user->J); CHKERRQ(ierr);
  
  /* solvers */
  ierr = SNESCreate(user->comm,&user->snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(user->snes,user->R,FormResidual,(void*)user);
  ierr = SNESSetJacobian(user->snes,user->J,user->J,FormJacobian,(void*)user);
  ierr = SNESSetFromOptions(user->snes);CHKERRQ(ierr);
  ierr = SNESSetUp(user->snes);

  /* timestep viewer */
  sprintf(filename,"%s_%s",user->param->filename,"ts.csv");
  ierr = PetscViewerASCIIOpen(user->comm,filename,&user->timestep_table);
  
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CleanUpDataStructures"
PetscErrorCode CleanUpDataStructures(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  /* display parameters again */
  ierr = PetscBagView(user->bag,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  /* clean up */
  ierr = VecDestroy(&user->X);CHKERRQ(ierr);
  ierr = VecDestroy(&user->Xo);CHKERRQ(ierr);
  ierr = VecDestroy(&user->R);CHKERRQ(ierr);
  ierr = MatDestroy(&user->J);CHKERRQ(ierr);
  ierr = SNESDestroy(&user->snes);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&user->timestep_table);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "DoSolve"
PetscErrorCode DoSolve(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  Parameter      *par = user->param;
  PetscInt       n_next_out = par->nout;
  PetscReal      lnR, F;
  SNESConvergedReason reason;
  PetscFunctionBeginUser;
  ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  while (par->n < par->ns && par->t < par->tmax && F < par->Fmax) {

    /* solve timestep */
    ierr = SNESSolve(user->snes,NULL,user->X);CHKERRQ(ierr);

    /* check solve outcome */
    ierr = SNESGetConvergedReason(user->snes,&reason);CHKERRQ(ierr);
    if (reason<0) { ierr = WriteOutput(user); CHKERRQ(ierr); break; }

    /* process timestep */
    par->t += par->dt; par->n++;
    ierr = VecCopy(user->X,user->Xo); CHKERRQ(ierr);
    ierr = TimestepTableAddEntry(user);CHKERRQ(ierr);
    ierr = VecGetValues(user->X,1,&par->ni,&lnR);CHKERRQ(ierr);
    F = 1 - pow(exp(lnR),3);
    ierr = PetscPrintf(user->comm,"Step: %d/%d, Time: %.5g/%.2g, F: %.4g/%.4g\n",par->n,par->ns,par->t,par->tmax,F,par->Fmax);CHKERRQ(ierr);
    if (par->n >= n_next_out) { ierr = WriteOutput(user); CHKERRQ(ierr); n_next_out += par->nout; }
    ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
    
  }
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "TimestepTableAddEntry"
PetscErrorCode TimestepTableAddEntry(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  Parameter      *par = user->param;
  PetscInt       Nv = 4+N_ODES, i, ind[Nv];
  PetscReal      val[Nv];
  PetscFunctionBeginUser;
  ierr = PetscViewerASCIIPrintf(user->timestep_table,"%d\t",par->n);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(user->timestep_table,"%7.7f\t",par->t);CHKERRQ(ierr);
  ind[0] = 0;  ind[1] = 1;
  for (i=2; i<Nv; i++) { ind[i] = par->ni + i - 4; }
  ierr = VecGetValues(user->X, Nv, ind, val);
  ierr = PetscViewerASCIIPrintf(user->timestep_table,"%7.7f\t",0.5*(val[0]+val[1]));CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(user->timestep_table,"%7.7f\t",0.5*(val[2]+val[3]));CHKERRQ(ierr);
  for (i=4; i<Nv; i++) { ierr = PetscViewerASCIIPrintf(user->timestep_table,"%7.7f\t",val[i]);CHKERRQ(ierr); }
  ierr = PetscViewerASCIIPrintf(user->timestep_table,"\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
