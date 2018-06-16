CFLAGS   = 
FFLAGS   = 
CPPFLAGS =
FPPFLAGS =
PROG     = superheat
OBJECTS  = ${PROG}.o 

all:	${PROG}
	${RM} ${OBJECTS}

debug:  ${PROG}

cleanup:
	-@${RM} ${OBJECTS}
	-@${RM} ${PROG}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

${PROG}: ${OBJECTS} chkopts
	 ${RM} ${PROG} ${CFLAGS}
	-${CLINKER} -o ${PROG} ${OBJECTS} ${PETSC_SNES_LIB} ${PETSC_CONTRIB}
