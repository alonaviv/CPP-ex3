CC = g++
FLAGS = -Wextra -Wall -Wvla -pthread -std=c++11 


all: Matrix.hpp.gch
	$(CC) $(FLAGS) Test.cpp Complex.cpp -o Test
	./Test

Matrix.hpp.gch: Matrix.hpp
	$(CC) $(FLAGS) Matrix.hpp
