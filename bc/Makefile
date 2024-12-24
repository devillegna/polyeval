

CC ?= clang
LD ?= clang


PROJ ?= avx2
#PROJ ?= ref
#PROJ ?= neon

OS := $(shell uname -s)
ARCH := $(shell uname -m)

ifneq (,$(filter $(ARCH),arm64 aarch64))
PROJ    = neon
endif

LIBNAME = fftbc


CFLAGS := -O3 -std=c11 -Wall -Wextra -Wpedantic -fno-omit-frame-pointer  #-Werror
LDFLAGS  := $(LDFLAGS)
INCPATH  = -I/usr/local/include -I/opt/local/include -I/usr/include -I./src
LIBS    = -lcrypto


EXT_SRC_DIRS  =

ifdef ASAN
CFLAGS  += -fsanitize=address -fno-sanitize-recover=all
LDFLAGS += -fsanitize=address -fno-sanitize-recover=all
endif


ifeq ($(PROJ),ref)

EXT_SRC_DIRS  += ./src/ref
INCPATH      += -I./src/ref

else ifeq ($(PROJ),avx2)

EXT_SRC_DIRS  += ./src/avx2
INCPATH      += -I./src/avx2
CFLAGS       += -mavx2 -mpclmul

else ifeq ($(PROJ),neon)

EXT_SRC_DIRS  += ./src/neon
INCPATH      += -I./src/neon

endif


TESTINCPATH  := -I/usr/local/include -I/opt/local/include -I/usr/include -I./src -I./benchmark -I./unit-tests 



#############################

CSRC        := $(wildcard $(EXT_SRC_DIRS)/*.c)
CSRC        += $(wildcard src/*.c)
SRC_O       := $(CSRC:.c=.o)
SRC_O_NODIR := $(notdir $(SRC_O))
OBJS   = $(SRC_O_NODIR)

#############################

OS := $(shell uname -s)
ARCH := $(shell uname -m)
ifeq  ($(OS), Darwin)
ifeq  ($(ARCH), arm64)
CFLAGS    +=  -D_APPLE_SILICON_
OBJS      += m1cycles.o
endif
endif


#############################

.INTERMEDIATE:  $(OBJS)

.PHONY: all clean


all: lib bc_32-test bc_512-test

%-test:  %-test.o lib
	$(LD) $(LDFLAGS) $(LIBPATH) -L. -o $@ $<  $(LIBS) -l$(LIBNAME)

%-benchmark:  %-benchmark.o lib
	$(LD) $(LDFLAGS) $(LIBPATH) -L. -o $@ $<  $(LIBS) -l$(LIBNAME)

%.o: unit-tests/%.c
	$(CC) $(CFLAGS) $(TESTINCPATH) -o $@ -c $<

%.o: benchmark/%.c
	$(CC) $(CFLAGS) $(TESTINCPATH) -o $@ -c $<

m1cycles.o: benchmark/m1cycles.c
	$(CC) $(CFLAGS) $(TESTINCPATH) -o $@ -c $<

lib: $(OBJS)
	ar rcs lib$(LIBNAME).a $^

%.o: src/%.c
	$(CC) $(CFLAGS) $(INCPATH) -o $@ -c $<

%.o: src/%.S
	$(CC) $(CFLAGS) $(INCPATH) -o $@ -c $<

define GEN_O
%.o: $(1)/%.c
	$(CC) $(CFLAGS) $(INCPATH) -o $$@ -c $$<
endef
$(foreach dir, $(EXT_SRC_DIRS), $(eval $(call GEN_O,$(dir))))


%.S: %.c
	$(CC) $(CFLAGS) -S -c  -o$@ $^


#############################

clean:
	-rm *.o *.s *.q *.a test speed *-test

