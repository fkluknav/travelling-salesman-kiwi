.PHONY: all clean valgrind test cppcheck

BINARY=run

all: test

$(BINARY): tsp.cc settings.h temperature.h ca_temp.h
	g++ -ffast-math -O4 -std=gnu++0x -Wall -Wextra -static -pthread -o $(BINARY) tsp.cc

test: $(BINARY)
	./$(BINARY) < data_15.txt
#time cat ../cut_data/data_70.txt | ./$(BINARY) | head --lines=2

clean:
	rm -f $(BINARY) $(BINARY)-valgrind  k.tsp k.sol Nk.*

#can not compile with -static
valgrind:
	g++ -ffast-math -ggdb3 -O4 -std=gnu++0x -Wall -Wextra -pthread -o $(BINARY)-valgrind tsp.cc
	cat ../3data | valgrind --trace-children=yes ./$(BINARY)-valgrind
	rm $(BINARY)-valgrind

cppcheck:
	cppcheck --enable=all ./tsp.cc
