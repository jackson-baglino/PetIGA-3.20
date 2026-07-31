#ifndef ENCELADUS_MAIN_H
#define ENCELADUS_MAIN_H

// Standard library includes
#include <math.h>
#include <mpi.h>

// PETSc and PetIGA includes
#include <petsc/private/tsimpl.h>
#include <petsc/private/snesimpl.h>
#include "petiga.h"

// Project-specific includes
#include "enceladus_types.h"
#include "material_properties.h"
#include "assembly.h"
#include "monitoring.h"
#include "initial_conditions.h"
#include "snes_convergence.h"

#endif // ENCELADUS_MAIN_H