CC = g++
CFLAGS = -Wall -g -std=c++11
DFLAGS = -lcrypto -pthread

all: brute_force

brute_force: brute_force.cpp
	$(CC) $(CFLAGS) $(DFLAGS) brute_force.cpp -o brute_forcer -L/ssl/lib -I/ssl/include

clean:
	rm *.o brute_forcer
