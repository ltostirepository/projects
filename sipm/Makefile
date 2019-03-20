CFLAGS = -O
CC = g++
SRC = src/hello.cpp
OBJ = $(SRC:.cpp = .o)

NumTest: $(OBJ)
    $(CC) $(CFLAGS) -o NumTest $(OBJ)
clean:
    rm -f core *.o
