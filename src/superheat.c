#include "superheat.h"

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
/*-----------------------------------------------------------------------*/
{
  AppCtx         user;  
  PetscErrorCode ierr;
  
  /* initialize PETSc */
  PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL);
  ierr = SetUpParameters(&user);CHKERRQ(ierr);
  ierr = SetUpDataStructures(&user);CHKERRQ(ierr);

  ierr = SNESSolve(user.snes,NULL,user.X);CHKERRQ(ierr);
  
  ierr = CleanUpDataStructures(&user);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormResidual"
PetscErrorCode FormResidual(SNES snes, Vec X, Vec R, void *ptr)
/* ------------------------------------------------------------------- */
{
  PetscErrorCode  ierr;
  AppCtx          *user = (AppCtx*)ptr;
  PetscReal const *x;
  PetscReal       *res;
  PetscInt        i;
  
  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x); CHKERRQ(ierr);
  ierr = VecGetArray(R,&res); CHKERRQ(ierr);

  for (i=0; i<user->param->dofs; i++) { res[i] = x[i] - i; }
  
  ierr = VecRestoreArray(R,&res); CHKERRQ(ierr);
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
  PetscFunctionBegin;

  /* Set up data structures */
  user->comm = PETSC_COMM_WORLD;
  ierr = PetscBagCreate(user->comm,sizeof(Parameter),&(user->bag)); CHKERRQ(ierr); 
  ierr = PetscBagGetData(user->bag,(void**)&user->param);CHKERRQ(ierr);
  ierr = PetscBagSetName(user->bag,"par","(parameters for crystal superheating problem)");CHKERRQ(ierr);
  par = user->param;

  /* Register parameters */
  ierr = PetscBagRegisterInt(user->bag,&par->ni,100,"ni","Number of grid points"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->dofs,par->ni+N_ODES,"dofs","<DO NOT SET> Total number of degrees of freedom"); CHKERRQ(ierr);
  ierr = PetscBagRegisterInt(user->bag,&par->ns,10000,"ns","Number of time steps"); CHKERRQ(ierr);

  /* Display parameters */
  ierr = PetscPrintf(user->comm,"--------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscBagView(user->bag,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(user->comm,"--------------------------------------\n");CHKERRQ(ierr);
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
  PetscFunctionBegin;

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

  /* SNES */
  ierr = SNESCreate(user->comm,&user->snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(user->snes,user->R,FormResidual,(void*)user);
  ierr = SNESSetFromOptions(user->snes);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CleanUpDataStructures"
PetscErrorCode CleanUpDataStructures(AppCtx *user)
/*-----------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecDestroy(&user->X);CHKERRQ(ierr);
  ierr = VecDestroy(&user->Xo);CHKERRQ(ierr);
  ierr = VecDestroy(&user->R);CHKERRQ(ierr);
  ierr = SNESDestroy(&user->snes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
