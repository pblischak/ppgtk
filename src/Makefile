# Makefile for compiling PPGtk with Boost and OpenMP

# Set this to `yes` if you have OpenMP.
OPENMP ?= yes

# Set this to where you have Boost installed
BOOST_LIB = /usr/local/lib
BOOST_INC = /usr/local/include

EXE = ppgtk
OBJ = DataParser.o MbRandom.o Frequency.o Phi.o Genotype.o ModelFreqs.o ModelDiseq.o main.o
TXT_FILES = frequencies.txt phi.txt

CXX = g++
CXXFLAGS = -I${BOOST_INC} -O3 -g -std=c++11
LDFLAGS = -lboost_program_options

ifeq ($(strip $(OPENMP)), yes)
	CXXFLAGS += -fopenmp
endif

.PHONY : clean install uninstall linux

$(EXE) : $(OBJ)
	@test -d ../bin || mkdir -p ../bin
	$(CXX) $(CXXFLAGS) -L${BOOST_LIB} -o $(EXE) $^ $(LDFLAGS)
	@mv $(EXE) ../bin
	@printf "\nTo install ppgtk system-wide run: sudo make install\n\n"

linux : $(OBJ)
	@test -d ../bin || mkdir -p ../bin
	$(CXX) $(CXXFLAGS) -L${BOOST_LIB} -Wl,-rpath=${BOOST_LIB} -o $(EXE) $^ $(LDFLAGS)
	@mv $(EXE) ../bin
	@printf "\nTo install ppgtk system-wide run: sudo make install\n\n"

$(OBJ) : %.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean :
	@printf "\n Removing object (*.o) files from src/ and executable files from bin/...\n\n"
	@rm -ir $(OBJ) $(TXT_FILES) ../bin

install :
	@printf "\nCopying ppgtk to /usr/local/bin...\n\n"
	@cp ../bin/$(EXE) /usr/local/bin
	@print "\nTo uninstall, type: sudo make uninstall\n\n"

uninstall :
	@printf "\nRemoving ppgtk from /usr/local/bin...\n\n"
	@rm -i /usr/local/bin/$(EXE)
