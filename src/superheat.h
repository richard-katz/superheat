#include "petsc.h"

#define FNAME_LENGTH  120
#define N_ODES 3

typedef struct {
  PetscInt  ni, ns, dofs;// number of grid points, max number of time steps
  PetscReal R0;          // initial radius
  PetscReal t, dt, tmax; // current time, time-step, maximum time
  char      filename[FNAME_LENGTH];
} Parameter;

typedef struct {
  Parameter     *param;
  PetscBag      bag;
  Vec           X, Xo, R;
  Mat           J;
  SNES          snes;
  MPI_Comm      comm;
} AppCtx;

PetscErrorCode SetUpParameters(AppCtx*);
PetscErrorCode SetUpDataStructures(AppCtx*);
PetscErrorCode CleanUpDataStructures(AppCtx*);
PetscErrorCode FormResidual(SNES, Vec, Vec, void*);
