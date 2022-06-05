# name of the executable to be generated
PROG = acotsp
# space-delimited list of header files
HDRS = acotsp.h util.h
# space-delimited list of source files
SRCS = acotsp.c main.c util.c

# directory to store object files
OBJDIR = object
# names of object files
OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(SRCS)) 

# name of the compiler
CC = gcc
# additional compiler flags to pass in
CFLAGS = -Wall --std=c99
# libraries for the linker
LIBS = -lm

####################
# Compiling rules: #
####################
# WARNING: *must* have a tab before each definition

# invoked when "make" is run
all : $(OBJDIR) $(PROG)

# links object files into executable
$(PROG) : $(OBJS)
	$(CC) $(CFLAGS) $^ -o $(PROG) $(LIBS)

# compiles source files into object files
object/%.o : %.c $(HDRS)
	$(CC) -c $(CFLAGS) $< -o $@ $(LIBS)

# creates directory to store object files
$(OBJDIR) :
	mkdir -p $@/

# cleans up object files and executable
# type "make clean" to use
# Note: you can add your own commands to remove other things (e.g. output files)
clean:
	rm -rf object/
	rm -f $(PROG)
