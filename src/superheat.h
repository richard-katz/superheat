#include "petsc.h"

#define FNAME_LENGTH  120
#define N_ODES 3

typedef struct {
  /* numerical model parameters */
  PetscInt  ni, dofs, N; // grid points for Cs, total degrees of freedom, output file number
  PetscInt  ns, n, nout; // max number of time steps, step number, output step interval
  PetscReal t, dt, tmax; // current time, time-step, maximum time
  PetscBool make_output;
  char      filename[FNAME_LENGTH];

  /* physical parameters */
  PetscReal K;      // Partition coefficient
  PetscReal decmpr; // Dimensionless decompression rate
  PetscReal St;     // Stefan number
  PetscReal phi;    // Dynamic liquid fraction (upper bound)
  PetscBool infinite_diffusion; // For case of homogeneous grain
} Parameter;

typedef struct {
  Parameter     *param;
  PetscBag      bag;
  Vec           X, Xo, R;
  Mat           J;
  SNES          snes;
  MPI_Comm      comm;
  PetscViewer   timestep_table;
} AppCtx;

PetscErrorCode SetUpParameters(AppCtx*);
PetscErrorCode SetUpDataStructures(AppCtx*);
PetscErrorCode SetUpInitialGuess(AppCtx*);
PetscErrorCode DoSolve(AppCtx*);
PetscErrorCode TimestepTableAddEntry(AppCtx*);
PetscErrorCode CleanUpDataStructures(AppCtx*);
PetscErrorCode WriteOutput(AppCtx*);
PetscErrorCode FormResidual(SNES, Vec, Vec, void*);
PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void*);
