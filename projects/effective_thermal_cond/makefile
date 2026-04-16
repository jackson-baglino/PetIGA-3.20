###############################################################################
# Makefile — effective_k_ice_homog
#
# Builds the effective thermal conductivity homogenization executable using
# PETSc and PetIGA.  Source code is split across several modules under src/.
#
# Targets:
#   make          — optimized build (-O3)
#   make debug    — debug build (-g3 -O0, ASan-friendly)
#   make clean    — remove executable and object files
#
# Required environment variables:
#   PETSC_DIR, PETSC_ARCH, PETIGA_DIR
###############################################################################

TARGETS = effective_k_ice_homog

all: $(TARGETS)
.PHONY: all

# ----------------------------------------------------------------------------
# Build flags
# ----------------------------------------------------------------------------
ifeq ($(BUILD),debug)
  CFLAGS = -g3 -O0 -Wall -Wextra -Wno-unused-parameter
else
  CFLAGS = -O3   -Wall -Wextra -Wno-unused-parameter
endif

debug:
	$(MAKE) BUILD=debug all
.PHONY: debug

# ----------------------------------------------------------------------------
# Compiler / linker
# ----------------------------------------------------------------------------
CC      = mpicc
LDFLAGS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc \
          -L$(PETIGA_DIR)/$(PETSC_ARCH)/lib -lpetiga \
          -lmpi -lmetis -lparmetis -lsuperlu_dist

# ----------------------------------------------------------------------------
# Include paths
# ----------------------------------------------------------------------------
INCLUDES = -I$(PETSC_DIR)/include \
           -I$(PETSC_DIR)/include/petsc/private \
           -I$(PETIGA_DIR)/include \
           -I./include
CFLAGS  += $(INCLUDES)

# ----------------------------------------------------------------------------
# PetIGA configuration (provides CLINKER, PETIGA_LIB, PETSC_LIB)
# ----------------------------------------------------------------------------
include $(PETIGA_DIR)/lib/petiga/conf/variables
include $(PETIGA_DIR)/lib/petiga/conf/rules
include $(PETIGA_DIR)/lib/petiga/conf/test

# ----------------------------------------------------------------------------
# Module object files
# ----------------------------------------------------------------------------
OBJS_HOMOG = obj/effective_k_ice_homog.o \
             obj/material.o              \
             obj/field_init.o            \
             obj/env_config.o            \
             obj/assembly.o              \
             obj/solver.o                \
             obj/io_thermal.o            \
             obj/setup.o

# ----------------------------------------------------------------------------
# Compilation rule
# ----------------------------------------------------------------------------
obj/%.o: src/%.c
	@mkdir -p obj
	$(CC) $(CFLAGS) -c $< -o $@

# ----------------------------------------------------------------------------
# Link rule
# ----------------------------------------------------------------------------
effective_k_ice_homog: $(OBJS_HOMOG)
	$(CLINKER) -o $@ $^ $(PETIGA_LIB) $(PETSC_LIB)

# ----------------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------------
clean::
	@echo "Cleaning..."
	@rm -f $(TARGETS)
	@find obj -type f -name '*.o' -delete 2>/dev/null || true
	@echo "Done."
