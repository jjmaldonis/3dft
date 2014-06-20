

FC = ifort
FCFLAGS = -O2 -openmp -fpp
#FCFLAGS = -g -debug all -check all -implicitnone -warn all -fpp

.PHONY: ft ift stdev clean
.SUFFIXES: .f90 .o
# Set default ways of creating .o files from .f90 and .c files
%.o: %.f90
	${FC} -c $< ${FCFLAGS}

ft: APP = 3dft
ft: SRC = gfx.f90 model_v2.f90 scattering_factors.f90 3dft.f90
ft: OBJ = $(addsuffix .o, $(basename $(SRC)))
ft: ${OBJ}
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

ift: APP = inverse3dft
ift: SRC = gfx.f90 model_v2.f90 scattering_factors.f90 inverse3dft.f90
ift: OBJ = $(addsuffix .o, $(basename $(SRC)))
ift: ${OBJ}
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

stdev: APP = stdev
stdev: SRC = gfx.f90 model_v2.f90 scattering_factors.f90 stdev.f90
stdev: OBJ = $(addsuffix .o, $(basename $(SRC)))
stdev: ${OBJ}
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

clean:
	rm -f *.mod *.o ${APP}
