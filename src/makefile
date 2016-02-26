#ifort -O2 -openmp -fpp -c 2dft.f90 && ifort gfx.o model_v2.o scattering_factors.o 2dft.o -O2 -openmp -fpp -o 2dft

FC = ifort
FCFLAGS = -O2 -openmp -fpp
#FCFLAGS = -g -debug all -check all -implicitnone -warn all -fpp -openmp

SRC = gfx.f90 model_v2.f90 scattering_factors.f90
OBJ = $(addsuffix .o, $(basename $(SRC)))

all: ft ift stdev

.PHONY: ft ift stdev clean all
.SUFFIXES: .f90 .o
# Set default ways of creating .o files from .f90 and .c files
%.o: %.f90
	${FC} -c $< ${FCFLAGS}

ft: ${OBJ} ${SRC} 3dft.o
	${FC} ${OBJ} 3dft.o ${FCFLAGS} -o ../3dft

ift: ${OBJ} ${SRC} inverse3dft.o
	${FC} ${OBJ} inverse3dft.o ${FCFLAGS} -o ../inverse3dft

stdev: ${OBJ} ${SRC} stdev.o
	${FC} ${OBJ} stdev.o ${FCFLAGS} -o ../stdev

clean:
	rm -f *.mod *.o 3dft inverse3dft stdev
