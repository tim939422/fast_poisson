# Compiler and flags
FC = gfortran
FFLAGS = -O2
FFTW_LIB = /Users/duosifan/Codes/fftw/3.3.10/lib
FFTW_INC = /Users/duosifan/Codes/fftw/3.3.10/include
LIBS = -lfftw3
INCLUDES = -I$(FFTW_INC)

EXE = poisson_1d.exe

$(EXE): main.o fftw.o m_poisson_solver.o m_constants.o
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ $^ -L$(FFTW_LIB) $(LIBS)

m_constants.o: m_constants.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c m_constants.f90

fftw.o: fftw.f90 m_constants.o
	$(FC) $(FFLAGS) $(INCLUDES) -c fftw.f90

m_poisson_solver.o: m_poisson_solver.f90 m_constants.o fftw.o
	$(FC) $(FFLAGS) $(INCLUDES) -c m_poisson_solver.f90

main.o: main.f90 fftw.o m_poisson_solver.o m_constants.o
	$(FC) $(FFLAGS) $(INCLUDES) -c main.f90

# Clean rule: remove executables and object files
.PHONY: clean
clean:
	rm -f $(EXE) *.o