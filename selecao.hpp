/*
*  selecao.hpp
*
*  Módulo responsável pelas reotinas de seleção
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 08/07/2017
*  Datas de Modificações:
*
*/

#ifndef _SELECAO_H
#define	_SELECAO_H

#include "individuo.hpp"
#include "populacao.hpp"
#include "utils.hpp"
#include "tsp.hpp"
#include "arqlog.hpp"

using namespace std;

class TSelecao
{
	public:
		static const int tipoSorteio   = 0;
	   static const int tipoTorneioK  = 1;
		static const int tipoRoleta    = 2;
		static const int tipoRoletaInv = 3;
		
	private:
      TMapaGenes *VP_Mapa;
      TArqLog *VP_ArqSaida;
		int VP_TipoSelecao;  //0 - sorteio
		                     //1 - torneio
									//2 - roleta
									//3 - roleta inversa
	
	public:
      //Apenas leitura
      TMapaGenes *getMapa();
      TArqLog *getArqLog();

	   TSelecao (TMapaGenes *mapa, TArqLog *arqSaida, int tipoSelecao);
		~TSelecao ();

		int getTipo();
		void setTipo(int val);

      unsigned processa(TPopulacao *populacao, unsigned exceto);
		unsigned processa(TPopulacao *populacao, unsigned exceto, int k);
		unsigned processa(TPopulacao *populacao);
		
	private:
      unsigned torneioK  (TPopulacao *populacao, int k);
		unsigned sorteio   (TPopulacao *populacao);
		unsigned roleta    (TPopulacao *populacao);
		unsigned roletaInv (TPopulacao *populacao);
};

#endif