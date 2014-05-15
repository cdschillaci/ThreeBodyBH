CC = g++

OPTIMIZE= -DNDEBUG -O3

CFLAGS = -c -fopenmp -I../eigen/ -Wall -g $(OPTIMIZE) 
LFLAGS = -lm -lpthread -fopenmp -Wall -g

OBJECTS = main.o gamma.o basis.o HOstate.o energies.o energies2.o energies3.o heff.o heff3.o

all: main 

main: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main

main.o: main.cpp basis.h energies3.h 
	$(CC) $(CFLAGS) main.cpp 

gamma.o: gamma.h gamma.cpp
	$(CC) $(CFLAGS) gamma.cpp 

basis.o: basis.h basis.cpp HOstate.h
	$(CC) $(CFLAGS) basis.cpp

HOstate.o: HOstate.h HOstate.cpp
	$(CC) $(CFLAGS) HOstate.cpp

heff.o: basis.h gamma.h heff.h heff.cpp
	$(CC) $(CFLAGS) heff.cpp

heff3.o: basis.h gamma.h heff.h heff3.h heff3.cpp
	$(CC) $(CFLAGS) heff3.cpp

energies.o: heff.h energies.h energies.cpp
	$(CC) $(CFLAGS) energies.cpp

energies2.o: heff.h energies.h energies2.h energies2.cpp
	$(CC) $(CFLAGS) energies2.cpp

energies3.o: heff.h energies.h energies3.h energies3.cpp
	$(CC) $(CFLAGS) energies3.cpp

clean: 
	rm *.o main
