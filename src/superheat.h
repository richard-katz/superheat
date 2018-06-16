#include "petsc.h"

#define FNAME_LENGTH  120

typedef struct {
  PetscInt  ni, ns;      // number of grid points, max number of time steps
  PetscReal R0;          // initial radius
  PetscReal t, dt, tmax; // current time, time-step, maximum time
  char      filename[FNAME_LENGTH];
} Parameter;

typedef struct {
  Parameter     *param;
  PetscBag      bag;
  Vec           X, Xo, R;
  SNES          snes;
  MPI_Comm      comm;
} AppCtx;
