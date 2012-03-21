
CXX=CC -O3 -fastsse -static

GSL_BASE=$(GSL_DIR)
GSL_INCL=-I$(GSL_BASE)/include
GSL_LDFLAGS=-L$(GSL_BASE)/lib
GSL_LIBS=-lgsl

TARGET = make_model

OBJECTS = $(TARGET).o

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
