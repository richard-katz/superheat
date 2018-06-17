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
  ierr = SetUpDataStructures(&user);CHKERRQ(ierr);

  /* form IC, solve and output */
  ierr = SetUpInitialGuess(&user);CHKERRQ(ierr);
  ierr = DoSolve(&user);CHKERRQ(ierr);
  ierr = WriteOutput(&user);CHKERRQ(ierr);

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
  PetscReal const *x, *xo;
  PetscReal       dt = user->param->dt;
  PetscReal       dr, r, R, *res, delsq, Rdot, Cdot, mvcrd;
  PetscInt        i, is, ie, iR;
  
  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(X,&x); CHKERRQ(ierr);
  ierr = VecGetArrayRead(user->Xo,&xo); CHKERRQ(ierr);
  ierr = VecGetArray(Res,&res); CHKERRQ(ierr);
  dr = 1./(user->param->ni-2);

  /* radius ODE */
  iR = user->param->dofs - N_ODES;
  R = x[iR];  Rdot  = (x[iR] - xo[iR])/dt;
  res[iR] = Rdot + 1; // FIX THIS
  
  /* diffusion boundary conditions */
  is = 0; ie = user->param->ni-1;
  res[is] = x[is] - x[is+1];     is++;
  res[ie] = x[ie] + x[ie-1] - 0; ie--;

  /* diffusion PDE */
  for (i=is; i<=ie; i++) {
    r = dr*(i-0.5);
    Cdot  = (x[i] - xo[i])/dt;
    delsq = ((x[i-1] - 2*x[i] + x[i+1])/dr/dr + (x[i+1] - x[i-1])/r/dr)/R/R;
    mvcrd = r*Rdot*(x[i+1] - x[i-1])/(2*dr);
    res[i] = mvcrd + delsq - Cdot;
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
  PetscReal const *x, *xo;
  PetscInt        row, col[4], is, ie, rowR;
  PetscReal       A[4], dr, r, dt, R, Rdot, delsq;
  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(X,&x); CHKERRQ(ierr);
  ierr = VecGetArrayRead(user->Xo,&xo); CHKERRQ(ierr);

  is = 0; ie = user->param->ni-1;
  dr = 1./(user->param->ni-2);
  dt = user->param->dt;

  /* radius ODE */
  rowR = user->param->dofs - N_ODES;
  R = x[rowR];  Rdot = (x[rowR] - xo[rowR])/dt; 
  col[0] = rowR; A[0] = 1/dt;
  ierr = MatSetValues(J,1,&rowR,1,col,A,INSERT_VALUES);CHKERRQ(ierr);

  /* diffusion boundary conditions */
  col[0]=is;   A[0] = +1;
  col[1]=is+1; A[1] = -1;
  ierr = MatSetValues(J,1,&is,2,col,A,INSERT_VALUES);CHKERRQ(ierr); is++; 
  col[0]=ie-1; A[0] = +1;
  col[1]=ie;   A[1] = +1;
  ierr = MatSetValues(J,1,&ie,2,col,A,INSERT_VALUES);CHKERRQ(ierr); ie--;
  
  /* diffusion PDE */
  for (row=is; row<=ie; row++) {
    r = dr*(row-0.5);
    col[0] = row-1; col[1] = row; col[2] = row+1; col[3] = rowR;
    /*left*/ A[0] = (1/dr - 1/r)/dr/R/R - r*Rdot/(2*dr);
    /*cntr*/ A[1] = -2/dr/dr/R/R - 1/dt;
    /*rght*/ A[2] = (1/dr + 1/r)/dr/R/R + r*Rdot/(2*dr);
    delsq = ((x[row-1] - 2*x[row] + x[row+1])/dr/dr + (x[row+1] - x[row-1])/r/dr)/R/R;    
    /*R   */ A[3] = r*(x[row+1] - x[row-1])/(2*dr)/dt - 2*delsq/R;
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

  /* Register parameters */
  ierr = PetscBagRegisterInt(user->bag,&par->ni,100,"ni","Number of grid points"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->dofs,par->ni+N_ODES,"dofs","<DO NOT SET> Total number of degrees of freedom"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->ns,100,"ns","Number of time steps"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->t,0,"t","Time");CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(user->bag,&par->dt,1e-3,"dt","Time-step size");CHKERRQ(ierr);
  ierr = PetscBagRegisterString(user->bag,&par->filename,FNAME_LENGTH,"test","filename","Name of output file");CHKERRQ(ierr);

  /* Display parameters */
  ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscBagView(user->bag,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "WriteOutput"
PetscErrorCode WriteOutput(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscFunctionBeginUser;
  ierr = PetscPrintf(user->comm," Generating output file: \"%s\"\n",user->param->filename);
  ierr = PetscViewerBinaryOpen(user->comm,user->param->filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
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
  PetscFunctionBeginUser;
  ierr = VecSet(user->X,1);CHKERRQ(ierr);
  ierr = VecSet(user->Xo,1);CHKERRQ(ierr);
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
  ierr = MatSeqAIJSetPreallocation(user->J,4,NULL); CHKERRQ(ierr); // FIX THIS
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
  PetscInt       n;
  PetscFunctionBeginUser;
  for (n=0; n<user->param->ns; n++) {
    ierr = SNESSolve(user->snes,NULL,user->X);CHKERRQ(ierr);
    ierr = VecCopy(user->X,user->Xo); CHKERRQ(ierr);
    user->param->t += user->param->dt;
    ierr = PetscPrintf(user->comm,"-----------------------------------------\n");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
