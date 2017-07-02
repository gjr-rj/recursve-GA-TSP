/*
*  mutacao.hpp
*
*  Módulo responsável pelo mutação de um indivíduo
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 14/05/2017
*  Datas de Modificações:
*
*/

#ifndef _MUTACAO_H
#define	_MUTACAO_H

#include <vector>
#include "individuo.hpp"
#include "utils.hpp"
#include "tsp.hpp"
#include "arqlog.hpp"

class TMutacao
{
   private:
      TMapaGenes *VP_Mapa;
      TArqLog *VP_ArqSaida;

      int VP_mutacao;
   public:

      //Apenas leitura
      TMapaGenes *getMapa();
      TArqLog *getArqLog();


      TMutacao (TMapaGenes *mapa, TArqLog *arqSaida, int tipoMutacao);
      int processa (TIndividuo *individuo);

   private:
      //Mutações
      int EM(TIndividuo *individuo);
      int SM(TIndividuo *individuo);
      int twoOPT(TIndividuo *individuo);
      int twoOPT2(TIndividuo *individuo);
      int SIM(TIndividuo *individuo);
      int threeOPT(TIndividuo *individuo);
      int threeOPT2(TIndividuo *individuo);
	  
	  //Métodos auxiliares
	  double calcValReverso (TGene *g1, TGene *g2);
};

#endif
