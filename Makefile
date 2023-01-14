BUILD_DIR = gen
SOURCE_DIR = src

OMPCC = gcc -fopenmp
MPICC = mpicc
NVCC = nvcc

CC_FLAGS = -O3
CC_FLAGS += -Wall -Wextra
NVCC_FLAGS = -O3
LIBS = -lm

ifeq ($(DEBUG), 1)
CC_FLAGS += -DDEBUG
endif

all: $(BUILD_DIR)/prime $(BUILD_DIR)/feynman $(BUILD_DIR)/moldyn \
$(BUILD_DIR)/prime-mpi $(BUILD_DIR)/feynman-mpi $(BUILD_DIR)/moldyn-mpi $(BUILD_DIR)/moldyn-mw-mpi \
$(BUILD_DIR)/prime-cuda

$(BUILD_DIR)/prime: $(SOURCE_DIR)/dz1z1.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/feynman: $(SOURCE_DIR)/dz1z3.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/moldyn: $(SOURCE_DIR)/dz1z5.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(OMPCC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/prime-mpi: $(SOURCE_DIR)/dz2z1.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(MPICC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/feynman-mpi: $(SOURCE_DIR)/dz2z2.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(MPICC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/moldyn-mpi: $(SOURCE_DIR)/dz2z3.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(MPICC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/moldyn-mw-mpi: $(SOURCE_DIR)/dz2z4.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(MPICC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/prime-cuda: $(SOURCE_DIR)/dz4z1.cu $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(NVCC) $(NVCC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
