CC= g++
CFLAGS=  -Wall -std=c++0x -g -I/usr/local/include -c
CFLAGS1= -L/usr/local/lib
OBJCS = main.o F2.o Data.o Steer.o Blad.o
NAME= F2.exe
all: main

main: $(OBJCS)
	$(CC) $(CFLAGS1) $(OBJCS) -lgsl -lgslcblas -lm -o $(NAME)
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
Blad.o: Blad.cpp
	$(CC) $(CFLAGS) Blad.cpp
Data.o: Data.cpp
	$(CC) $(CFLAGS) Data.cpp
Steer.o: Steer.cpp
	$(CC) $(CFLAGS) Steer.cpp
F2.o: F2.cpp
	$(CC) $(CFLAGS) F2.cpp
clean:
	rm *.o
