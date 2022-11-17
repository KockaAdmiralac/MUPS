BUILD_DIR = gen
SOURCE_DIR = src

CC = gcc

CC_FLAGS = -fopenmp -O3
CC_FLAGS += -Wall -Wextra
LIBS = -lm

ifeq ($(DEBUG), 1)
CC_FLAGS += -DDEBUG
endif

all: $(BUILD_DIR)/prime $(BUILD_DIR)/feynman $(BUILD_DIR)/moldyn

$(BUILD_DIR)/prime: $(SOURCE_DIR)/dz1z1.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(CC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/feynman: $(SOURCE_DIR)/dz1z3.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(CC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR)/moldyn: $(SOURCE_DIR)/dz1z5.c $(SOURCE_DIR)/util.c | $(BUILD_DIR)
	$(CC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
