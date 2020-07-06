CONSERVATIVE_FLAGS = -std=c++11 -Wall -Wextra -pedantic -O3
DEBUGGING_FLAGS = -g
CFLAGS = $(CONSERVATIVE_FLAGS) $(DEBUGGING_FLAGS)

mosaicism: main.o DiscriminantFunctions.o FilterFunctions.o
	g++ -o mosaicism main.o DiscriminantFunctions.o FilterFunctions.o

FilterFunctions.o: FilterFunctions.h FilterFunctions.cpp
	g++ $(CFLAGS) -c FilterFunctions.cpp

DiscriminantFunctions.o: DiscriminantFunctions.h DiscriminantFunctions.cpp
	g++ $(CFLAGS) -c DiscriminantFunctions.cpp

main.o: main.cpp DiscriminantFunctions.h FilterFunctions.h
	g++ $(CFLAGS) -c main.cpp

.PHONY: clean
clean:
	rm -f *.o mosaicism?  *~ mosaicism
