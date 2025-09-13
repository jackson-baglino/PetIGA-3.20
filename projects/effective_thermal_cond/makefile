###############################################################################
# Makefile for effective_k_ice and effective_k_ice_full
#
# This makefile compiles two executables using PETSc and PetIGA:
#   1. effective_k_ice       – modular build with separate source files
#   2. effective_k_ice_full  – single-file build for easier debugging
#
# It includes configuration files from PetIGA which set additional variables
# (CLINKER, PETIGA_LIB, PETSC_LIB) required for proper linking.
#
# Build Modes:
#   - Optimized Build (default):
#       make all         – Compiles with high optimization flags (e.g., -O3)
#   - Debug Build:
#       make debug       – Compiles with debugging flags (e.g., -g3 -O0)
#
# Usage:
#   make               – Builds both executables in optimized mode
#   make debug         – Builds both executables in debug mode
#   make clean         – Removes executables and object files
#
# Environment variables required:
#   PETSC_DIR, PETSC_ARCH, PETIGA_DIR
###############################################################################

#
# ----------------------------------------------------------------------------
# List of executables to build
# ----------------------------------------------------------------------------
TARGETS = effective_k_ice_homog
# 					effective_k_ice_full \
					effective_k_ice_homog

all: $(TARGETS)
.PHONY: all

# ----------------------------------------------------------------------------
# Choose build flags based on BUILD variable (debug vs optimized)
# ----------------------------------------------------------------------------
ifeq ($(BUILD),debug)
  CFLAGS = -g3 -O0 -Wall -Wextra -Wno-unused-parameter $(INCLUDES)
  FFLAGS = -g3 -O0 -Wall -Wextra -fcheck=all
else
  CFLAGS = -O3   -Wall -Wextra -Wno-unused-parameter $(INCLUDES)
  FFLAGS = -O3   -Wall -Wextra
endif

debug:
	$(MAKE) BUILD=debug all
.PHONY: debug

# ----------------------------------------------------------------------------
# Compiler Settings
# ----------------------------------------------------------------------------
# CC: MPI C compiler wrapper
CC       = mpicc
# ----------------------------------------------------------------------------
# Linker Settings
# ----------------------------------------------------------------------------
# LDFLAGS: Link against PETSc, PetIGA, MPI, and other dependencies
LDFLAGS  = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc \
           -L$(PETIGA_DIR)/$(PETSC_ARCH)/lib -lpetiga \
           -lmpi -lmetis -lparmetis -lsuperlu_dist

# ----------------------------------------------------------------------------
# Include Paths
# ----------------------------------------------------------------------------
# Header search paths for PETSc, PetIGA, and local headers
INCLUDES = -I$(PETSC_DIR)/include \
           -I$(PETSC_DIR)/include/petsc/private \
           -I$(PETIGA_DIR)/include \
           -I./include
CFLAGS  += $(INCLUDES)

# ----------------------------------------------------------------------------
# PetIGA Configuration
# ----------------------------------------------------------------------------
# Include PetIGA-provided variables, rules, and test settings
include $(PETIGA_DIR)/lib/petiga/conf/variables
include $(PETIGA_DIR)/lib/petiga/conf/rules
include $(PETIGA_DIR)/lib/petiga/conf/test

# ----------------------------------------------------------------------------
# Source Files
# ----------------------------------------------------------------------------
# Gather all .c source files under src/
SRCS         := $(wildcard src/*.c)
# ----------------------------------------------------------------------------
# Common Modules
# ----------------------------------------------------------------------------
# Exclude the two main programs from the common object list
SRCS_COMMON  := $(filter-out src/effective_k_ice.c  \
                             src/effective_k_ice_full.c \
														 src/effective_k_ice_homog.c, $(SRCS))
OBJS_COMMON  := $(patsubst src/%.c, obj/%.o, $(SRCS_COMMON))

# ----------------------------------------------------------------------------
# Main Program Objects
# ----------------------------------------------------------------------------
# Object files for each main program
OBJ_MAIN_ICE := obj/effective_k_ice.o
OBJ_MAIN_FULL := obj/effective_k_ice_full.o
OBJ_MAIN_HOMOG := obj/effective_k_ice_homog.o

OBJS_ICE     := $(OBJ_MAIN_ICE)   $(OBJS_COMMON)
OBJS_FULL    := $(OBJ_MAIN_FULL)  # standalone, contains all code already

# ----------------------------------------------------------------------------
# Compilation Rule
# ----------------------------------------------------------------------------
# Compile each .c in src/ into an object in obj/
obj/%.o: src/%.c
	@mkdir -p obj
	$(CC) $(CFLAGS) -c $< -o $@

# ----------------------------------------------------------------------------
# Link Rule
# ----------------------------------------------------------------------------
# Link the corresponding object files to produce the executable
effective_k_ice:   $(OBJS_ICE)
	$(CLINKER) -o $@ $^ $(PETIGA_LIB) $(PETSC_LIB)

effective_k_ice_full: $(OBJS_FULL)
	$(CLINKER) -o $@ $^ $(PETIGA_LIB) $(PETSC_LIB)

effective_k_ice_homog: $(OBJ_MAIN_HOMOG)
	$(CLINKER) -o $@ $^ $(PETIGA_LIB) $(PETSC_LIB)

# ----------------------------------------------------------------------------
# Clean Target
# ----------------------------------------------------------------------------
# Remove executables and all object files
clean::
	@echo "Cleaning up..."
	@rm -f $(TARGETS)
	@find obj -type f -name '*.o' -delete
	@echo "Clean"

# ----------------------------------------------------------------------------
# Debug Build Target
# ----------------------------------------------------------------------------
# Invoke make with BUILD=debug to produce debug versions (no optimization)