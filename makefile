#ifort -O2 -openmp -fpp -c 2dft.f90 && ifort gfx.o model_v2.o scattering_factors.o 2dft.o -O2 -openmp -fpp -o 2dft

FC = ifort
FCFLAGS = -O2 -openmp -fpp
#FCFLAGS = -g -debug all -check all -implicitnone -warn all -fpp

SRC = gfx.f90 model_v2.f90 scattering_factors.f90
OBJ = $(addsuffix .o, $(basename $(SRC)))
ft: APP = 3dft
ft: SRC = gfx.f90 model_v2.f90 scattering_factors.f90 3dft.f90
ft: OBJ = $(addsuffix .o, $(basename $(SRC)))
ift: APP = inverse3dft
ift: SRC = gfx.f90 model_v2.f90 scattering_factors.f90 inverse3dft.f90
ift: OBJ = $(addsuffix .o, $(basename $(SRC)))
stdev: APP = stdev
stdev: SRC = gfx.f90 model_v2.f90 scattering_factors.f90 stdev.f90
stdev: OBJ = $(addsuffix .o, $(basename $(SRC)))


.PHONY: ft ift stdev clean all
.SUFFIXES: .f90 .o
# Set default ways of creating .o files from .f90 and .c files
%.o: %.f90
	${FC} -c $< ${FCFLAGS}

ft: ${OBJ} ${SRC} 3dft.o
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

ift: ${OBJ} ${SRC} inverse3dft.o
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

stdev: ${OBJ} ${SRC} stdev.o
	${FC} ${OBJ} ${FCFLAGS} -o ${APP}

clean:
	rm -f *.mod *.o 3dft inverse3dft stdev
