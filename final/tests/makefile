#We try to detect the OS we are running on, and adjust commands as needed
CLEANUP = rm -f
MKDIR = mkdir -p
TARGET_EXTENSION=.out

#Tool Definitions
CC=gcc
CFLAGS=-I. -I$(PATHU) -DTEST
CFLAGS+=`gsl-config --cflags --libs`

#Path Definitions
PATHU = ../../tools/unity/
PATHS =
PATHT =
PATHB = build/

#Files We Are To Work With
SRC = test_of_unity.c unity.c
OBJ = $(patsubst %.c,$(PATHB)%.o,$(SRC))
DEP = $(PATHU)unity.h $(PATHU)unity_internals.h
TGT = $(PATHB)test$(TARGET_EXTENSION)

test: $(PATHB) $(TGT)
	./$(TGT)

$(PATHB)%.o:: $(PATHS)%.c $(DEP)
	$(CC) -c $(CFLAGS) $< -o $@

$(PATHB)%.o:: $(PATHU)%.c $(DEP)
	$(CC) -c $(CFLAGS) $< -o $@

$(TGT): $(OBJ)
	gcc -o $@ $^

clean:
	$(CLEANUP) $(PATHB)*.o
	$(CLEANUP) $(TGT)

$(PATHB):
	$(MKDIR) $(PATHB)

all: clean test

.PHONY: all
.PHONY: clean
.PHONY: test
