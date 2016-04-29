# Makefile for compiling PPGtk with Boost and OpenMP
# Modified from the Makefile for EEMS (https://github.com/dipetkov/eems)
# and SVDquartets (http://www.stat.osu.edu/~lkubatko/software/SVDquartets/).

# Set this to `no` if you don't have OpenMP.
OPENMP ?= yes

# Set this to where you have Boost installed
BOOST_LIB = /usr/local/lib
BOOST_INC = /usr/local/include

EXE = ppgtk
OBJ = CmdLineParser.o DataParser.o MbRandom.o main.o

CXX = g++
CXXFLAGS = -I${BOOST_INC} -O3 -DNDEBUG -std=c++11
LDFLAGS = \
	-lboost_system \
	-lboost_program_options \
	-lboost_filesystem

ifeq ($(strip $(OPENMP)), yes)
	CXXFLAGS += -fopenmp
endif

.PHONY : all linux darwin clean install uninstall docs

all :
	@echo "Use 'make linux' on a Linux machine, or 'make darwin' on a Mac."

linux : $(OBJ)
	$(CXX) $(CXXFLAGS) -L${BOOST_LIB} -Wl,-R${BOOST_LIB} $(LDFLAGS) -o $(EXE) $^
	@mv $(EXE) ../bin
	@printf "\nTo install ppgtk system-wide run: sudo make install\n\n"

darwin : $(OBJ)
	@mkdir ../bin
	$(CXX) $(CXXFLAGS) -L${BOOST_LIB} $(LDFLAGS) -o $(EXE) $^
	@mv $(EXE) ../bin
	@printf "\nTo install ppgtk system-wide run: sudo make install\n\n"

$(OBJ) : %.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean :
	@printf "\n Removing object (*.o) files from src/ and executable files from bin/...\n\n"
	@rm -ir $(OBJ) ../bin

install :
	@printf "\nCopying ppgtk to /usr/local/bin...\n\n"
	@cp ../bin/ppgtk /usr/local/bin
	@print "\nTo uninstall, type: sudo make uninstall\n\n"

uninstall :
	@printf "\nRemoving ppgtk from /usr/local/bin...\n\n"
	@rm -i /usr/local/bin/ppgtk
