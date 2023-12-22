## CS130 Fall 2013: Computer Graphics
## Lab 1 Makefile
##
## This file does not need to be modified
#########################################

LIBS = -lglut -lGLU -lGL -lXmu -lXext -lXi -lX11 -O3 -lm
CC = g++

## Global header files
INCLUDE =

## Object files and executables
MAIN_OUT = render

## Requirements for each command
MAIN_REQS = main.cpp

## Targets to compile for each command
MAIN_TARGETS = main.cpp

all: main

## Main 
main: $(MAIN_REQS)
	$(CC) -ggdb $(MAIN_TARGETS) $(LIBS) -o $(MAIN_OUT)

clean:
	rm -f *~ render *.out
