CC = gcc
CFLAGS = -std=c11 -O3 -pedantic -Wall -fopenmp
LDLIBS = -lm -llapacke -lgsl -lgslcblas -fopenmp

SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin

SRCS := fock.c gpe_radial.c quick_function.c rho.c state.c util.c testing.c examples.c
SRCS := $(addprefix $(SRC_DIR)/,$(SRCS))
OBJS := $(SRCS:$(SRC_DIR)/%=$(OBJ_DIR)/%.o)

all: $(OBJ_DIR) run

run: $(OBJS) $(OBJ_DIR)/main.c.o
	$(CC) $^ $(LDLIBS) -o run

$(OBJ_DIR)/main.c.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.c.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm obj/*.o
