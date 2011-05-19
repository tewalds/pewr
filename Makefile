.PHONY: clean fresh run gendeps profile


LDFLAGS   += -lfftw3 -lm
OBJECTS		= pewr.o

ifdef DEBUG
	CPPFLAGS	+= -g3 -Wall
else
	LDFLAGS     += -fopenmp
	CPPFLAGS    += -fopenmp
	CPPFLAGS    += -march=native
	CPPFLAGS	+= -O3 -funroll-loops -ffast-math -ftree-vectorize -ftree-loop-im -Wall
#	CPPFLAGS	+= -O3 -fopenmp -funroll-loops -ffast-math -funsafe-math-optimizations -ftree-vectorize -ftree-loop-im -mfpmath=sse -mmmx -msse -msse2 -msse3 -m3dnow
endif

#profile with callgrind, works well with DEBUG mode
ifdef PROFILE
	CFLAGS		+= -pg
endif

#For profile directed optimizations. To use:
# 1. compile with PROFILE_GEN
# 2. run the program, generates .gcno and .gcda
# 3. compile with PROFILE_USE
ifdef PROFILE_GEN
	CFLAGS		+= -fprofile-generate
endif
ifdef PROFILE_USE
	CFLAGS		+= -fprofile-use
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


