#include "petsc.h"

#define FNAME_LENGTH  120
#define N_ODES 1

typedef struct {
  PetscInt  ni, ns, dofs;// number of grid points, max number of time steps
  PetscReal t, dt, tmax; // current time, time-step, maximum time
  PetscBool make_output;
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
PetscErrorCode SetUpInitialGuess(AppCtx*);
PetscErrorCode DoSolve(AppCtx*);
PetscErrorCode CleanUpDataStructures(AppCtx*);
PetscErrorCode WriteOutput(AppCtx*);
PetscErrorCode FormResidual(SNES, Vec, Vec, void*);
PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void*);
