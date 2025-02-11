

CC = clang
LD = clang

LIB3NAME = fftbc
LIB3PATH = ./bc

LIB1NAME = gf264btfy
LIB1PATH = ./btfy/gf264

LIB2NAME = gf232btfy
LIB2PATH = ./btfy/gf232

LIB0NAME = gf264dencoder
LIB0PATH = ./dencoder/gf264


CFLAGS   := -O3 -std=c11 -Wall -Wextra -Wpedantic -fno-omit-frame-pointer  #-Werror
INCPATH  := -I/usr/local/include -I/opt/local/include -I/usr/include -I./src
LDFLAGS  := $(LDFLAGS) -L$(LIB0PATH) -L$(LIB1PATH)  -L$(LIB2PATH) -L$(LIB3PATH)
LIBS     := -lcrypto -l$(LIB0NAME) -l$(LIB1NAME)  -l$(LIB2NAME) -l$(LIB3NAME)


PROJ ?= avx2
#PROJ ?= ref
#PROJ ?= neon

OS := $(shell uname -s)
ARCH := $(shell uname -m)

ifneq (,$(filter $(ARCH),arm64 aarch64))
PROJ    = neon
endif

EXT_SRC_DIRS  =

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
ifeq  ($(ARCH), aarch64)
CFLAGS    +=  -flax-vector-conversions -march=armv8-a+crypto -mfpu=neon
else
CFLAGS    +=  -flax-vector-conversions -mcpu=apple-m1
endif
endif

ifdef ASAN
CFLAGS  += -fsanitize=address -fno-sanitize-recover=all
LDFLAGS += -fsanitize=address -fno-sanitize-recover=all
endif


TESTINCPATH  := -I/usr/local/include -I/opt/local/include -I/usr/include -I./include -I./src -I./benchmark -I./unit-tests -I$(LIB0PATH)/src -I$(LIB1PATH)/src  -I$(LIB2PATH)/src -I$(LIB3PATH)/src

INCPATH += $(TESTINCPATH)

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
#OBJS      += m1cycles.o
endif
endif


#############################

.INTERMEDIATE:  $(OBJS)

.PHONY: all clean


all: lib gf264_eval-test gf232_eval-test

lib: $(LIB0PATH)/lib$(LIB0NAME).a $(LIB1PATH)/lib$(LIB1NAME).a  $(LIB2PATH)/lib$(LIB2NAME).a $(LIB3PATH)/lib$(LIB3NAME).a

$(LIB0PATH)/lib$(LIB0NAME).a:
	cd $(LIB0PATH) && $(MAKE) PROJ=$(PROJ) lib

$(LIB1PATH)/lib$(LIB1NAME).a:
	cd $(LIB1PATH) && $(MAKE) PROJ=$(PROJ) lib

$(LIB2PATH)/lib$(LIB2NAME).a:
	cd $(LIB2PATH) && $(MAKE) PROJ=$(PROJ) lib

$(LIB3PATH)/lib$(LIB3NAME).a:
	cd $(LIB3PATH) && $(MAKE) PROJ=$(PROJ) lib

%-test: $(OBJS) %-test.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^  $(LIBS)

%-benchmark: $(OBJS) %-benchmark.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^  $(LIBS)

%.o: unit-tests/%.c
	$(CC) $(CFLAGS) $(TESTINCPATH) -o $@ -c $<

%.o: benchmark/%.c
	$(CC) $(CFLAGS) $(TESTINCPATH) -o $@ -c $<

m1cycles.o: benchmark/m1cycles.c
	$(CC) $(CFLAGS) $(TESTINCPATH) -o $@ -c $<

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

