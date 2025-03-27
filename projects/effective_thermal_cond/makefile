# =============================
# 🔹 Project Configuration
# =============================

# Define the main target (executable name)
TARGETS = effective_k_ice

# Define the directory structure
SRC_DIR = src         # Directory for source files (.c)
OBJ_DIR = obj         # Directory for compiled object files (.o)
INC_DIR = include     # Directory for header files (.h)
SCRIPT_DIR = scripts  # Directory for scripts (e.g., shell scripts)

# =============================
# 🔹 Compiler & Linker Configuration
# =============================

# Compiler flags
CFLAGS = -g3 -O0 -Wall -Wextra -Wno-unused-parameter \
         -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/petsc/private \
         -I$(PETIGA_DIR)/include -I$(INC_DIR)

# Include paths (modularized for easy modification)
INCLUDES = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/petsc/private \
           -I$(PETIGA_DIR)/include -I$(INC_DIR)
CFLAGS += $(INCLUDES)

# Explanation of CFLAGS:
# - `-g3`     : Enables debugging information.
# - `-O0`     : Disables optimization (useful for debugging).
# - `-Wall`   : Enables most warnings.
# - `-Wextra` : Enables extra warnings.
# - `-Wno-unused-parameter` : Suppresses warnings about unused function parameters.
# - `-I...`   : Specifies directories to search for header files.

# Linker flags - specify library paths and libraries to link
LDFLAGS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -L$(PETIGA_DIR)/lib -lpetiga

# Explanation of LDFLAGS:
# - `-L<dir>` : Specifies directories where libraries are located.
# - `-lpetsc` : Links the PETSc library.
# - `-lpetiga` : Links the PetIGA library.

# Include configuration files from PetIGA (these set up additional variables and rules)
include $(PETIGA_DIR)/lib/petiga/conf/variables
include $(PETIGA_DIR)/lib/petiga/conf/rules

# =============================
# 🔹 Source & Object Files
# =============================

# Find all .c files in the source directory
SRCS := $(wildcard $(SRC_DIR)/*.c)

# Replace the source file paths with corresponding object file paths
OBJS := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# =============================
# 🔹 Build Rules
# =============================

# Ensure object directory exists before building
all: $(OBJ_DIR) $(TARGETS)

# $(OBJ_DIR):
# 	@mkdir -p $(OBJ_DIR)

# Compile each .c file into an .o file
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Explanation:
# - `$<` : The input file (source file).
# - `$@` : The output file (object file).
# - `$(CC)` : Uses the default C compiler.

# Link all object files to create the final executable
$(TARGETS): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# Explanation:
# - `$@` : The output file (the final executable).
# - `$^` : All prerequisites (object files).

# =============================
# 🔹 Utility Commands
# =============================

.PHONY: all clean run scripts

# Clean up compiled files
clean:
	-@$(RM) $(TARGETS)
	-@$(RM) $(OBJ_DIR)/*.o
	@rmdir $(OBJ_DIR) 2>/dev/null || true

# Explanation:
# - `.PHONY: clean` : Ensures `make clean` always runs, even if there's a file named "clean".
# - `$(RM)` : Removes files (typically `rm -f`).
# - `$(OBJ_DIR)/*.o` : Deletes all compiled object files.
# - `@rmdir $(OBJ_DIR) 2>/dev/null || true` : Silently remove `obj/` if empty.

# Run the simulation (calls a shell script)
run: all
	bash $(SCRIPT_DIR)/run_simulation.sh

# Explanation:
# - Ensures the code is compiled before running the simulation.
# - Executes a shell script that sets up and runs the simulation.

# Execute all scripts in the scripts directory
scripts:
	@for script in $(wildcard $(SCRIPT_DIR)/*.sh); do \
		bash $$script; \
	done

# Explanation:
# - `$(wildcard $(SCRIPT_DIR)/*.sh)` : Finds all `.sh` scripts.
# - `bash $$script` : Runs each script found.
# - The `@` suppresses echoing of the command itself.
