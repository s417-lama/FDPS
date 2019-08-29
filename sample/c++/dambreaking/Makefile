PS_PATH = -I ../../../src/
MYTH_PATH = $(HOME)/opt/massivethreads

CC = g++
#CC = mpicxx
CFLAGS = -O3 -ffast-math -funroll-loops -std=c++11 -march=native
# CFLAGS = -O0 -g -std=c++11

# # Enable OpenMP
# CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

# Enable MassiveThreads
CFLAGS += -DPARTICLE_SIMULATOR_TASK_PARALLEL -I$(MYTH_PATH)/include
LIBS = -lmyth -L $(MYTH_PATH)/lib -Wl,-R$(MYTH_PATH)/lib

#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = sph.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CC) -c $< $(CFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

clean:
	rm -f *.out *.o

distclean: clean
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
