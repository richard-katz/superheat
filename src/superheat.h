#include "petsc.h"

#define FNAME_LENGTH  120
#define N_ODES 2

typedef struct {
  /* numerical model parameters */
  PetscInt  ni, ns, dofs;// number of grid points, max number of time steps
  PetscReal t, dt, tmax; // current time, time-step, maximum time
  PetscBool make_output;
  char      filename[FNAME_LENGTH];

  /* physical parameters */
  PetscReal K;     // Partition coefficient
  PetscReal Pdot;  // Dimensionless decompression rate
  PetscReal St;    // Stefan number
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
