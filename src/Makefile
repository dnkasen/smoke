
#CXX=CC -O3 -fastsse -static
CXX=mpicxx -O3 -openmp
GSL_DIR=/home5/baker/friesen/gsl-1.15/intel-12.1.2

GSL_BASE=$(GSL_DIR)
GSL_INCL=-I$(GSL_BASE)/include
GSL_LDFLAGS=-L$(GSL_BASE)/lib
GSL_LIBS=-lgsl

TARGET = smoke

OBJECTS = main.o grid.o transport.o spectrum.o radioactive.o

$(TARGET): $(OBJECTS)
	$(CXX) $(GSL_LDFLAGS) -o $(TARGET) $(OBJECTS) $(GSL_LIBS)

.cc.o: 
	$(CXX) $(GSL_INCL) -c -o $@ $<

.c.o: 
	cc -I. -c -o $@ $<

clean:
	rm -f $(OBJECTS)

realclean: clean
	rm -f $(TARGET)
