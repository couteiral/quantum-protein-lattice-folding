#
#  PAULID ANNEALER MAKEFILE
#


CXX       = mpic++
CXXFLAGS  = -O3 --std=c++11 $(PETSC_FC_INCLUDES)
OBJECTS   = main.o matrix_builder.o options.o hamiltonian.o timer.o

include  $(PETSC_DIR)/lib/petsc/conf/variables
include  $(PETSC_DIR)/lib/petsc/conf/rules


annealer:   $(OBJECTS) chkopts
	@echo "========= PauliD annealer compilation ========="
	@echo "CXX    = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "========= PauliD annealer compilation ========="
	-${CLINKER} -o annealer ${OBJECTS} ${PETSC_LIB}
	${RM} $(OBJECTS)
