#
#  PAULID EIGENSOLVER MAKEFILE
#


CXX       = mpic++
CXXFLAGS  = -O3 --std=c++11 $(PETSC_FC_INCLUDES)
OBJECTS   = main.o matrix_builder.o options.o hamiltonian.o timer.o

include $(PETSC_DIR)/lib/petsc/conf/variables
#include $(PETSC_DIR)/lib/petsc/conf/rules
include $(SLEPC_DIR)/lib/slepc/conf/slepc_variables
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common


eigensolver:   $(OBJECTS) chkopts
	@echo "========= PauliD eigensolver compilation ========="
	@echo "CXX    = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "========= PauliD eigensolver compilation ========="
	-${CLINKER} -o eigensolver ${OBJECTS} ${PETSC_LIB} ${SLEPC_SYS_LIB}
	${RM} $(OBJECTS)


