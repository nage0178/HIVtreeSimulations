CC = gcc
vpath %.c src testing
vpath %.h src
CFLAGS=`pkg-config --cflags glib-2.0` -g -Wall -O3
LDADD=`pkg-config --libs glib-2.0`
GSLFLAGS=-lgsl -lgslcblas -lm

dna_sim: dna_sim.o string_utilities.o dna_sim_functions.o
	$(CC) $(CFLAGS) dna_sim.o dna_sim_functions.o string_utilities.o -o dna_sim $(LDADD) $(GSLFLAGS)

dna_sim.o: dna_sim.c dna_sim.h
	$(CC) $(CFLAGS) -c $<  $(GSLFLAGS)

dna_sim_functions.o: dna_sim_functions.c dna_sim.h
	$(CC) $(CFLAGS) -c $<

string_utilities.o: string_utilities.c string_utilities.h
	$(CC) $(CFLAGS) -c $<

test_dna_functions: dna_sim_functions.o test_dna_functions.o string_utilities.o dna_sim_functions.o
	$(CC) $(CFLAGS) test_dna_functions.o dna_sim_functions.o string_utilities.o -o test_dna_functions $(LDADD) $(GSLFLAGS)

test_dna_functions.o: test_dna_functions.c dna_sim.h
	$(CC) $(CFLAGS) -c $< $(GSLFLAGS)

compare_tree: compare_tree.o dna_sim_functions.o string_utilities.o
	$(CC) $(CFLAGS) compare_tree.o string_utilities.o dna_sim_functions.o -o compare_tree $(LDADD) $(GSLFLAGS)

compare_tree.o: compare_tree.c dna_sim.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) sasprintf.o
	$(RM) sasprintf
	$(RM) dna_sim
	$(RM) dna_sim.o
	$(RM) dna_sim_functions.o
	$(RM) string_utilities.o
	$(RM) test_dna_functions
	$(RM) compare_tree
	$(RM) compare_tree.o
