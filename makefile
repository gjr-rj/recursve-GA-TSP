LIB = lib/
EXE = tsp
CPP = g++
#CPPFLAGS =  -g -ggdb -Wall -Wextra -Wno-sign-compare -fopenmp -lm -std=c++11
#CPPFLAGS =  -Wall -Wextra -Wno-sign-compare -O2 -pipe -march=native -fopenmp -lm -std=c++11
CPPFLAGS =  -Wall `xml2-config --cflags --libs` -std=c++11
OBJ = main.o utils.o tipos.o tsp.o config.o ag.o arqlog.o populacao.o individuo.o mutacao.o cruzamento.o tabconversao.o selecao.o
all: 	$(OBJ)
	$(CPP) $(OBJ) -o  $(EXE) $(CPPFLAGS)
clean:
	$(RM) $(EXE)
	$(RM) $(OBJ)
main.o: main.cpp
utils.o: utils.cpp utils.hpp
tipos.o: tipos.cpp
tsp.o: tsp.cpp tsp.hpp
config.o: config.cpp config.hpp
ag.o: ag.cpp ag.hpp
arqlog.o: arqlog.cpp arqlog.hpp
populacao.o: populacao.cpp populacao.hpp
individuo.o: individuo.cpp individuo.hpp
mutacao.o: mutacao.cpp mutacao.hpp
cruzamento.o: cruzamento.cpp cruzamento.hpp
tabconversao.o: tabconversao.cpp tabconversao.hpp
selecao.o: selecao.cpp selecao.hpp