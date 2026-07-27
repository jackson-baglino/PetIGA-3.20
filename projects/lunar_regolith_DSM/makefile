TARGETS = permafrost

all: $(TARGETS)
.PHONY: all clean debug

###############################################################################
# Include PetIGA/PETSc configuration
# Provides: PCC, CLINKER, PETIGA_LIB, PETSC_LIB, CCPPFLAGS, etc.
###############################################################################
include $(PETIGA_DIR)/lib/petiga/conf/variables
include $(PETIGA_DIR)/lib/petiga/conf/rules
include $(PETIGA_DIR)/lib/petiga/conf/test

###############################################################################
# Include paths
# petsc/private is required for tsimpl.h / snesimpl.h used in NASA_main.h
###############################################################################
INCLUDES = -I$(PETSC_DIR)/include \
           -I$(PETSC_DIR)/include/petsc/private \
           -I$(PETIGA_DIR)/include \
           -I./include

###############################################################################
# Build mode
#   make           — optimized (default)
#   make debug     — full debug symbols, no optimisation
#   make BUILD=debug <target>  — same as above
###############################################################################
ifeq ($(BUILD),debug)
  OPT_FLAGS  = -g3 -O0
  WARN_FLAGS = -Wall -Wextra -Wno-unused-parameter
else
  # -march=native / -mtune=native: enable all M2 NEON/SIMD instructions
  # -funroll-loops: benefits tight Gauss-point loops in IGA assembly
  # -fomit-frame-pointer: frees one callee-save register in hot paths
  # -pipe: use pipes instead of temp files (faster compilation)
  OPT_FLAGS  = -O3 -march=native -mtune=native -funroll-loops \
               -fomit-frame-pointer -pipe
  WARN_FLAGS = -Wall -Wextra -Wno-unused-parameter
endif

# -MMD -MP: emit a .d dependency file alongside each .o so that editing
#           any header triggers the correct recompiles without make clean.
CFLAGS = $(OPT_FLAGS) $(WARN_FLAGS) $(INCLUDES) -MMD -MP

###############################################################################
# Source and object files
###############################################################################
SRCS := $(wildcard src/*.c)
OBJS := $(patsubst src/%.c, obj/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

###############################################################################
# Build rules
###############################################################################
obj/%.o: src/%.c
	@mkdir -p obj
	$(PCC) $(CFLAGS) -c $< -o $@

$(TARGETS): $(OBJS)
	$(CLINKER) -o $@ $^ $(PETIGA_LIB)

# Pull in auto-generated header dependencies (silently skipped if absent)
-include $(DEPS)

###############################################################################
# Clean
###############################################################################
clean::
	@echo "Cleaning..."
	@rm -f $(TARGETS)
	@rm -rf obj
	@echo "Done."

###############################################################################
# Debug build shorthand
###############################################################################
debug:
	$(MAKE) BUILD=debug all
