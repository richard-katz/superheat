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
  PetscReal       dr, r, R, *res, Rdot, Cdot, Vdot;
  PetscReal       CdotR, Cl, Cldot, Vl, Cs;
  PetscInt        i, iR, iCl, iV, iCs;
  
  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(X,&x); CHKERRQ(ierr);
  ierr = VecGetArrayRead(user->Xo,&xo); CHKERRQ(ierr);
  ierr = VecGetArray(Res,&res); CHKERRQ(ierr);
  dr = 1./(user->param->ni-2);

  /* extract ODE variables */
  iCs = user->param->ni - 1;            Cs = 0.5*(x[iCs] + x[iCs-1]);
  iR  = user->param->dofs - N_ODES;     R  = exp(x[iR]);  
  iCl = user->param->dofs - N_ODES + 1; Cl = x[iCl]; 
  iV  = user->param->dofs - N_ODES + 2; Vl = x[iV]; 

  /* radius ODE (incomplete) */
  Rdot  = (x[iR] - xo[iR])/dt;
  CdotR = (Cs - 0.5*(xo[iCs] + xo[iCs-1]))/dt;
  res[iR] = Rdot - (-par->decmpr - CdotR)/(par->St - (x[iCs] - x[iCs-1])/dr);

  /* liquid volume ODE (incomplete) */
  Vdot = (x[iV] - xo[iV])/dt;
  res[iV] = Vdot + 3*pow(R,3)*Rdot;
  
  /* liquid concentration ODE (incomplete) */
  Cldot = (x[iCl] - xo[iCl])/dt;
  res[iCl] = Vl*Cldot/pow(R,3)/3
           + Rdot*(par->K-1)*Cl
           + (x[iCs] - x[iCs-1])/dr/R/R;
  
  /* diffusion boundary conditions (complete) */
  res[0]   = x[0]   - x[1];
  res[iCs] = Cs - par->K*(Cl-1);

  /* diffusion PDE (complete) */
  for (i=1; i<iCs; i++) {
    r = dr*(i-0.5);  Cdot = (x[i] - xo[i])/dt;
    res[i] = -Cdot + r*Rdot*(x[i+1] - x[i-1])/(2*dr) +
           + ((x[i-1] - 2*x[i] + x[i+1])/dr/dr + (x[i+1] - x[i-1])/r/dr)/R/R;
  }
  
  ierr = VecRestoreArray(Res,&res); CHKERRQ(ierr);
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
  PetscReal       A[5], dr, r, dt, R, Rdot, delsq, Cl, Cldot, Vdot, Vl;
  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(X,&x); CHKERRQ(ierr);
  ierr = VecGetArrayRead(user->Xo,&xo); CHKERRQ(ierr);

  is = 0; ie = user->param->ni-1;
  dr = 1./(user->param->ni-2);
  dt = user->param->dt;

  /* extract ODE variables */
  iCs = user->param->ni - 1;
  iR  = user->param->dofs - N_ODES;     R  = exp(x[iR]);  
  iCl = user->param->dofs - N_ODES + 1; Cl = x[iCl]; 
  iV  = user->param->dofs - N_ODES + 2; Vl = x[iV]; 

  /* radius ODE (complete) */
  R = exp(x[iR]);  Rdot = (x[iR] - xo[iR])/dt; 
  /*R */       col[0] = iR;    A[0] = 1/dt;
  /*Cs left*/  col[1] = ie-1;  A[1] = +dr*(xo[iCs] + xo[iCs-1] - 2*(x[iCs]   + par->decmpr*dt) + par->St*dr)
				    / (2*dt*pow(par->St*dr - x[iCs] + x[iCs-1],2));
  /*Cs right*/ col[2] = ie;    A[2] = -dr*(xo[iCs] + xo[iCs-1] - 2*(x[iCs-1] + par->decmpr*dt) - par->St*dr)
				    / (2*dt*pow(par->St*dr - x[iCs] + x[iCs-1],2));
  ierr = MatSetValues(J,1,&iR,3,col,A,INSERT_VALUES);CHKERRQ(ierr);

  /* liquid volume ODE (incomplete) */
  Vl = x[iV];  Vdot = (x[iV] - xo[iV])/dt;
  /* Vl */      col[0] = iV;   A[0] = 1/dt;
  /* R  */      col[1] = iR;   A[1] = 3*pow(R,3)*(3*x[iR]*Rdot + 1/dt);
  ierr = MatSetValues(J,1,&iV,2,col,A,INSERT_VALUES);CHKERRQ(ierr);
  
  /* liquid concentration ODE (complete) */
  Cl = x[iCl]; Cldot = (x[iCl] - xo[iCl])/dt;
  /*Cl */      col[0] = iCl;  A[0] = Vl/dt/pow(R,3)/3 + Rdot*(par->K-1); 
  /*R  */      col[1] = iR;   A[1] = - Vl*Cldot/pow(R,3)*x[iR]
				     + (par->K-1)*Cl/dt
				     - 2*(x[ie]-x[ie-1])/dr/R/R*x[iR];
  /*Vl */      col[2] = iV;   A[2] = Cldot/pow(R,3)/3;
  /*Cs left*/  col[3] = ie-1; A[3] = -1/dr/R/R;
  /*Cs right*/ col[4] = ie;   A[4] = +1/dr/R/R;
  ierr = MatSetValues(J,1,&iCl,5,col,A,INSERT_VALUES);CHKERRQ(ierr);

  /* diffusion boundary conditions (complete) */
  col[0]=is;   A[0] = +1;
  col[1]=is+1; A[1] = -1;
  ierr = MatSetValues(J,1,&is,2,col,A,INSERT_VALUES);CHKERRQ(ierr);
  col[0]=ie-1; A[0] = +0.5;
  col[1]=ie;   A[1] = +0.5;
  col[2]=iCl;  A[2] = -par->K;
  ierr = MatSetValues(J,1,&ie,2,col,A,INSERT_VALUES);CHKERRQ(ierr);
  
  /* diffusion PDE (complete) */
  for (row=is+1; row<ie; row++) {
    r = dr*(row-0.5);
    col[0] = row-1; col[1] = row; col[2] = row+1; col[3] = iR;
    /*left*/ A[0] = (1/dr - 1/r)/dr/R/R - r*Rdot/(2*dr);
    /*cntr*/ A[1] = -2/dr/dr/R/R - 1/dt;
    /*rght*/ A[2] = (1/dr + 1/r)/dr/R/R + r*Rdot/(2*dr);
    delsq = ((x[row-1] - 2*x[row] + x[row+1])/dr + (x[row+1] - x[row-1])/r)/dr/R/R;    
    /*R   */ A[3] = r*(x[row+1] - x[row-1])/(2*dr)/dt - 2*delsq*x[iR];
    ierr = MatSetValues(J,1,&row,4,col,A,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
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
  ierr = PetscBagRegisterInt(user->bag,&par->ni,100,"ni","Number of grid points"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->dofs,par->ni+N_ODES,"dofs","<DO NOT SET> Total number of degrees of freedom"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->ns,200,"ns","Number of time steps"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->nout,10,"nout","Output step interval"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->n,0,"n","<DO NOT SET> Current time-step number"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->N,0,"N","<DO NOT SET> Current output frame number"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->t,0,"t","Time");CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->tmax,10,"tmax","Maximum time");CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->dt,1e-4,"dt","Time-step size");CHKERRQ(ierr);
  ierr = PetscBagRegisterString(user->bag,&par->filename,FNAME_LENGTH,"test","filename","Name of output file");CHKERRQ(ierr);

  /* Register physical parameters */
  ierr = PetscBagRegisterReal(user->bag,&par->K,1e-2,"K","Parition coefficient");CHKERRQ(ierr);  
  ierr = PetscBagRegisterReal(user->bag,&par->decmpr,1e-6,"decmpr","Dimensionless decompression rate");CHKERRQ(ierr);  
  ierr = PetscBagRegisterReal(user->bag,&par->St,0.1,"St","Stefan number");CHKERRQ(ierr);  
  
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
  sprintf(filename,"%s_%4.4d",user->param->filename,user->param->N++);
  ierr = PetscPrintf(user->comm,"Generating output file: \"%s\"\n",filename);
  ierr = PetscViewerBinaryOpen(user->comm,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = PetscBagView(user->bag,viewer);CHKERRQ(ierr); 
  ierr = VecView(user->X,viewer);CHKERRQ(ierr);
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
  ierr = VecSet(user->X,0);CHKERRQ(ierr);
  /* ln R */ ii[0] = user->param->dofs - N_ODES; vals[0] = 0;
  /* Cl   */ ii[1] = ii[0]+1;                    vals[1] = 1;
  /* Vl   */ ii[2] = ii[0]+2;                    vals[2] = 0.01;
  ierr = VecSetValues(user->X,N_ODES,ii,vals,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecCopy(user->X,user->Xo);CHKERRQ(ierr);
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
  ierr = VecDestroy(&user->X);CHKERRQ(ierr);
  ierr = VecDestroy(&user->Xo);CHKERRQ(ierr);
  ierr = VecDestroy(&user->R);CHKERRQ(ierr);
  ierr = MatDestroy(&user->J);CHKERRQ(ierr);
  ierr = SNESDestroy(&user->snes);CHKERRQ(ierr);
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
  PetscInt       n_next_out=0;
  SNESConvergedReason reason;
  PetscFunctionBeginUser;
  ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  while (par->n < par->ns && par->t < par->tmax) {
    ierr = SNESSolve(user->snes,NULL,user->X);CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(user->snes,&reason);CHKERRQ(ierr);
    if (reason<0) break;
    ierr = VecCopy(user->X,user->Xo); CHKERRQ(ierr);
    par->t += par->dt; par->n++;
    ierr = PetscPrintf(user->comm,"Step: %d/%d, Time: %.5g/%.2g\n",par->n,par->ns,par->t,par->tmax);CHKERRQ(ierr);
    if (par->n >= n_next_out) { ierr = WriteOutput(user); CHKERRQ(ierr); n_next_out += par->nout; }
    ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
