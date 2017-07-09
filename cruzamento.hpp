/*
*  cruzamento.hpp
*
*  Módulo responsável pelas reotinas de cruzamento
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 19/05/2017
*  Datas de Modificações:
*
*/

#ifndef _CRUZAMENTO_H
#define	_CRUZAMENTO_H

#include <set>
#include "selecao.hpp"
#include "individuo.hpp"
#include "populacao.hpp"
#include "utils.hpp"
#include "tsp.hpp"
#include "arqlog.hpp"

using namespace std;

class TCruzamento
{
	template<template <typename> class P = std::less > struct compare_pair_second 
	{
		template<class T1, class T2> bool operator()(const std::pair<T1, T2>& left, const std::pair<T1, T2>& right) 
		{
			return P<T2>()(left.second, right.second);
		}
	};

   private:
      TMapaGenes *VP_Mapa;
      TArqLog *VP_ArqSaida;

      int VP_cruzamento;
   public:

      //Apenas leitura
      TMapaGenes *getMapa();
      TArqLog *getArqLog();

      TCruzamento (TMapaGenes *mapa, TArqLog *arqSaida, int tipoCruzamento);
      vector <TIndividuo *>processa (TIndividuo *parceiro1, TIndividuo *parceiro2, TPopulacao *populacao, int maxGeracao, int geracao);

   private:
      vector <TIndividuo *>GSTX (TIndividuo *parceiro1, TIndividuo *parceiro2);
      vector <TIndividuo *>PMX  (TIndividuo *parceiro1, TIndividuo *parceiro2);
		vector <TIndividuo *>OX1  (TIndividuo *parceiro1, TIndividuo *parceiro2);
      vector <TIndividuo *>OX2  (TIndividuo *parceiro1, TIndividuo *parceiro2);
		vector <TIndividuo *>MOX  (TIndividuo *parceiro1, TIndividuo *parceiro2);
		vector <TIndividuo *>POS  (TIndividuo *parceiro1, TIndividuo *parceiro2);
      vector <TIndividuo *>CX   (TIndividuo *parceiro1, TIndividuo *parceiro2);
      vector <TIndividuo *>DPX  (TIndividuo *parceiro1, TIndividuo *parceiro2);
		vector <TIndividuo *>AP   (TIndividuo *parceiro1, TIndividuo *parceiro2);
		vector <TIndividuo *>MPX  (TIndividuo *parceiro1, TIndividuo *parceiro2);
		vector <TIndividuo *>HX   (TIndividuo *parceiro1, TIndividuo *parceiro2);
		vector <TIndividuo *>IO   (TIndividuo *parceiro1, TPopulacao *populacao);
		vector <TIndividuo *>MIO  (TIndividuo *parceiro1, TPopulacao *populacao, int maxGeracao, int geracao);
		vector <TIndividuo *>VR   (TIndividuo *parceiro1, TIndividuo *parceiro2, TPopulacao *populacao);
		vector <TIndividuo *>ER   (TIndividuo *parceiro1, TIndividuo *parceiro2);
};
#endif
