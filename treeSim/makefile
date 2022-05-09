CC = gcc
vpath %.c src testing
vpath %.h src
CFLAGS=-Wall -g -O3
LDADD=`pkg-config --libs glib-2.0`

HIV_treeSim: HIV_treeSim.o HIV_treeSim_functions.o HIV_treeSim.h
	$(CC) $(CFLAGS) HIV_treeSim.o HIV_treeSim_functions.o -o HIV_treeSim -lgsl -lgslcblas -lm

HIV_treeSim.o: HIV_treeSim.c HIV_treeSim.h
	$(CC) $(CFLAGS) -c $<

HIV_treeSim_functions.o: HIV_treeSim_functions.c  HIV_treeSim.h
	$(CC) $(CFLAGS) -c $< -std=gnu11

test_functions: test_functions.o HIV_treeSim_functions.o
		$(CC) $(CFLAGS) test_functions.o HIV_treeSim_functions.o -o test_functions -std=gnu11 -lgsl -lgslcblas -lm

test_functions.o: test_functions.c HIV_treeSim.h
		$(CC) $(CFLAGS) -c $< -std=gnu11

clean:
	$(RM) HIV_treeSim.o
	$(RM) HIV_treeSim
	$(RM) HIV_treeSim_functions.o
	$(RM) test_functions
	$(RM) test_functions.o
