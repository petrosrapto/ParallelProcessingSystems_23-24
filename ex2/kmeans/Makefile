.KEEP_STATE:

CC = gcc
CFLAGS = -Wall -Wextra -Wno-unused -O3
OMPFLAGS = -fopenmp $(CFLAGS)
# the -fopenmp instructs the gcc compiler to enable openmp directives in the \
source code, allowing the compiler to recognize and correctly compile code \
intented for parallel execution
LDFLAGS = -fopenmp
# when linking the final executable, the -fopenmp flag is essential to link \
OpenMP library with our program. The OpenMP runtime library contains code that \
runs the OpenMP directives at runtime, such as creating threads and handling parallel regions
SRC_DIR = src/
BUILD_DIR = build/
H_FILES = $(addprefix $(SRC_DIR), kmeans.h)
COMM_SRC = $(addprefix $(SRC_DIR), file_io.c util.c)

# Ensure the build directory exists
$(shell mkdir -p $(BUILD_DIR))

# Check if SRC_FILE variable is set. If not, print an error and exit
ifndef SRC_FILE
$(error No file given as an argument to make. As an example, submit to torque as follows: qsub -q parlab -v SRC_FILE=seq_kmeans run_on_queue.sh)
endif

# Conditional assignment of FLAGS
ifeq ($(SRC_FILE),seq_kmeans)
    FLAGS = $(CFLAGS)
else
    FLAGS = $(OMPFLAGS)
endif

all: $(BUILD_DIR)$(SRC_FILE)

$(BUILD_DIR)$(SRC_FILE): $(BUILD_DIR)main.o $(BUILD_DIR)file_io.o $(BUILD_DIR)util.o $(BUILD_DIR)$(SRC_FILE).o
	$(CC) $(FLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)main.o: $(SRC_DIR)main.c $(H_FILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)$(SRC_FILE).o: $(SRC_DIR)$(SRC_FILE).c $(COMM_SRC) $(H_FILES)
	$(CC) $(FLAGS) -c $< -o $@

$(BUILD_DIR)file_io.o: $(SRC_DIR)file_io.c
	$(CC) $(FLAGS) -c $< -o $@

$(BUILD_DIR)util.o: $(SRC_DIR)util.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)*

