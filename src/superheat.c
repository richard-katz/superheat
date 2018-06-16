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
  user.comm = PETSC_COMM_WORLD;
  ierr = PetscFinalize();
  return 0;
}
