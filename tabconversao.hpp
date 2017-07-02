/*
*  tabconversao.hpp
*
*  Módulo responsável pela conversão entre indivíduos.
          Auxilia a recursividade
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 21/05/2017
*  Datas de Modificações:
*
*/

#ifndef _TABECONVERSAO_H
#define	_TABECONVERSAO_H

#include <vector>
#include "individuo.hpp"
#include "populacao.hpp"
#include "tsp.hpp"

class TTabConversao
{
   private:
      vector<TTipoConversao> VP_tabConv;

   public:
      TTabConversao(TIndividuo *melhor, TIndividuo *select);
      ~TTabConversao();
	  int get_quantidade();
	  vector<TTipoConversao> get_vetConv();
	  void converte(TIndividuo *individuo, TIndividuo *indReduzido);
	  string toString();
};

#endif // _TABECONVERSAO_H
