


#APP = 3dft
#SRC = gfx.f90 model_v2.f90 scattering_factors.f90 3dft.f90

#APP = inverse3dft
#SRC = gfx.f90 model_v2.f90 scattering_factors.f90 inverse3dft.f90

APP = stdev
SRC = gfx.f90 model_v2.f90 scattering_factors.f90 stdev.f90

OBJ = $(addsuffix .o, $(basename $(SRC)))

FC = ifort
FCFLAGS = -O2 -openmp -fpp
#FCFLAGS = -g -debug all -check all -implicitnone -warn all -fpp

.SUFFIXES: .f90 .o
# Set default ways of creating .o files from .f90 and .c files
%.o: %.f90
	${FC} -c $< ${FCFLAGS}

default: ${OBJ}
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

clean:
	rm -f *.mod *.o ${APP}
