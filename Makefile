.PHONY: clean fresh run gendeps profile


LDFLAGS		+= -lm
OBJECTS		= pewr.o

ifdef USE_FLOATS
	CPPFLAGS	+= -DUSE_FLOATS
	LDFLAGS		+= -lfftw3f
else
	LDFLAGS		+= -lfftw3
endif

ifdef SINGLE_THREAD
	CPPFLAGS	+= -DSINGLE_THREAD
else
	LDFLAGS     += -fopenmp
	CPPFLAGS    += -fopenmp
endif

ifdef DEBUG
	CPPFLAGS	+= -g3 -Wall
else
	#CPPFLAGS    += -march=native
	CPPFLAGS	+= -O3 -funroll-loops -ffast-math -ftree-vectorize -ftree-loop-im -Wall
endif

#profile with callgrind, works well with DEBUG mode
ifdef PROFILE
	CPPFLAGS	+= -pg
endif

#For profile directed optimizations. To use:
# 1. compile with PROFILE_GEN
# 2. run the program, generates .gcno and .gcda
# 3. compile with PROFILE_USE
ifdef PROFILE_GEN
	CPPFLAGS	+= -fprofile-generate
endif
ifdef PROFILE_USE
	CPPFLAGS	+= -fprofile-use
endif


all: pewr

pewr: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

pewr.o: pewr.cpp Array.h time.h


gendeps:
	ls *.cpp -1 | xargs -L 1 cpp -M -MM

clean:
	rm -f pewr *.o

fresh: clean all

profile:
	valgrind --tool=callgrind ./pewr

