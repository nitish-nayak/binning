test: test.cpp
	g++ -std=c++17 -O3 -Wall -Wextra -pedantic test.cpp -o test
clean:
	rm test
