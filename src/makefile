################################################################################
#  qfeper Makefile                                                             #
#  Code author: Masoud Kazemi                                                  #
#  Latest update: 2017-06-09                                                   #
#  make [option] [COMP=compiler]                                               #
################################################################################


################################################################################
# Intel Fortran ifort ( tested in linux and osx)
################################################################################
ifeq ($(COMP),ifort)
        FC=             ifort
        FC_OPT=         -O3
        FFLAGS=         -ip -ipo ${FC_OPT}
        DEBUG=          -debug all -traceback -check all

################################################################################
# GNU gcc-gfortran
################################################################################
else ifeq ($(COMP),gcc)
        FC=             gfortran
        FC_OPT=         -O3
        FFLAGS=         -funroll-loops -ffree-line-length-none -fcray-pointer \
                        -static-libgfortran -fall-intrinsics -std=legacy \
                        -Wall -Wtabs -fstack-protector ${FC_OPT} -DG95=1

        FPP_FLG=        -cpp
        DEBUG=          -fbacktrace -g ${FC_OPT}

################################################################################
# If no COMP option is given assume gcc
################################################################################
else
        FC=             gfortran
        FC_OPT=         -O3
        FFLAGS=         -funroll-loops -ffree-line-length-none -fcray-pointer \
                        -static-libgfortran -fall-intrinsics -std=legacy \
                        -Wall -Wtabs -fstack-protector ${FC_OPT} -DG95=1

        FPP_FLG=        -cpp
        DEBUG=          -fbacktrace -g ${FC_OPT}
endif


################################################################################
# Display options
################################################################################
default: what

################################################################################
# Targets
################################################################################
all:	qfeper move1 move2

debug:
	@make FFLAGS="${DEBUG}" \
	qfeper move1 move2 

clean:
	-rm -f *.o  *.mod *.M *.kmo *.il *.oo

nuke:
	-rm -rf *.o  *.mod *.M *.kmo *.il *.oo qfeper ../bin ../obj 

qfeper: qfeper_pars.o qfeper_analyz.o qfeper.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@ 




################################################################################
# Object modules
################################################################################
qfeper_pars.o:qfeper_pars.f90
	${FC} ${FFLAGS} -c $<

qfeper_analyz.o:qfeper_analyz.f90
	${FC} ${FFLAGS} -c $<

qfeper.o:qfeper.f90 qfeper_pars.o qfeper_analyz.o
	${FC} ${FFLAGS} -c qfeper.f90

move1:
	mkdir -p ../obj ; mv *.o *.mod ../obj/

move2:
	mkdir -p ../bin ; mv qfeper ../bin/

################################################################################
# Build instructions
################################################################################
what:
	@echo "Use: make [target] [COMP=compiler]"
	@echo
	@echo "[target] is one of:"
	@echo "all      The qfeper program"
	@echo "clean    Clean all .o object files"
	@echo "nuke     Clean all .o files and also the qfeper executable"
	@echo
	@echo "Examples: make all COMP=ifort"
	@echo "          make all COMP=gcc"
	@echo "          make debug COMP=gcc"
################################################################################
