#- 
#- File    : Makefile
#- ------------------



# put your own GSL directory and compiler here if you want
GSL_DIR   = $(HOME)/gsl-1.15/intel-12.1.2
CXX   = mpicxx

#CXX   = CC


GSL_INC=$(GSL_DIR)/include
GSL_LIB=$(GSL_DIR)/lib

CCOPT = -I$(GSL_INC) -O3
CLOPT = -L$(GSL_LIB) -lgsl -lgslcblas


TARGET = smoke

OBJECTS = main.o grid.o transport.o spectrum.o radioactive.o

$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(CLOPT) $(FLIB)

.cc.o: 
	$(CXX) $(CCOPT) -c -o $@ $<

.c.o: 
	cc -I. -c -o $@ $<

clean:
	rm -f $(OBJECTS) smoke

