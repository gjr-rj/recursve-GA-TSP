/*
*  populacao.hpp
*
*  Módulo responsável tratamento da população
*  PPGI - NCE - UFRJ
*  Data Criação: 13/05/2017
*  Datas de Modificações:
*
*/

#ifndef _POPULACAO_H
#define	_POPULACAO_H

#include "tsp.hpp"
#include "arqlog.hpp"
#include "individuo.hpp"
#include <vector>
#include <algorithm>

class TPopulacao
{
   static int comparaMelhor(const void *Indi, const void *Indj);

   private:
      vector <TIndividuo *> VP_individuos;

      double VP_somaDistancias;
		double VP_somaDistanciasInv; // 1/VP_somaDistancias
      unsigned VP_tamanho;

      TMapaGenes *VP_Mapa;
      TArqLog *VP_ArqSaida;

   public:
      //Apenas leitura
      TMapaGenes *getMapa();
      TArqLog *getArqLog();
      unsigned get_tamanho();
      unsigned get_qtdeIndividuo();
      double get_soma_dist ();
		double get_soma_distInv ();
      void sub_dist_tot (double val);
      void soma_dist_tot (double val);

      //Métodos
      TPopulacao (int tamanho, TMapaGenes *mapa, TArqLog *arqSaida);
      ~TPopulacao ();

      void povoa();
      void povoa(vector<TTipoConversao> tabConv, int extraPadrao);

      TIndividuo *get_melhor();
      TIndividuo *get_pior();
      double distanciaMedia ();


      TIndividuo *get_individuo(unsigned indice);
		void troca (int indice1, int indice2);
      void add_individuo (TIndividuo *individuo);
		void add_individuo(vector <TIndividuo *>individuos);
      string toString ();
      void ordena();
};

#endif
