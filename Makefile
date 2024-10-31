all: seesaw_main

seesaw_main: seesaw_main.cpp seesaw.cpp seesaw.h
	g++ -I /opt/homebrew/include/eigen3 -I . -std=c++17 -lgsl -lgslcblas -Wall -o seesaw_main  seesaw_main.cpp seesaw.cpp
