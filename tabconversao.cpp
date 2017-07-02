/*
*  tabconversao.cpp
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

#include "tabconversao.hpp"
#include <iostream>

TTabConversao::TTabConversao(TIndividuo *melhor, TIndividuo *select)
{
   TGene *g;
   TGene *tmp;
   TTipoConversao item;

   bool primeiro = true;

   for(int i=0; i<melhor->get_qtdeGenes();i++)
   {
      g = melhor->get_por_indice(i);
      tmp = select->get_por_id(g->id);

      if (primeiro)
      {
         item.pr = g->id;
         item.prOrig = g->ori;
         item.prDest = g->dest;
         primeiro = false;
      }

      if ((g->prox->id != tmp->prox->id)||(g->prox->id==0))
      {
         item.ul = g->id;
         item.ulOrig = g->ori;
         item.ulDest = g->dest;
         VP_tabConv.push_back(item);
         primeiro = true;
      }
   }
}

TTabConversao::~TTabConversao()
{
   VP_tabConv.clear();
}

int TTabConversao::get_quantidade()
{
	return VP_tabConv.size();
}

vector<TTipoConversao> TTabConversao::get_vetConv()
{
   return VP_tabConv;
}

string TTabConversao::toString()
{
   string resultado = "";
   for (vector<TTipoConversao>::iterator i=VP_tabConv.begin(); i!=VP_tabConv.end(); ++i)
   {
      resultado += to_string((*i).pr);
	  resultado += "(";
	  resultado += to_string((*i).prOrig);
	  resultado += ",";
	  resultado += to_string((*i).prDest);
	  resultado += ") / ";
	  resultado += to_string((*i).ul);
	  resultado += "(";
	  resultado += to_string((*i).ulOrig);
	  resultado += ",";
	  resultado += to_string((*i).ulDest);
	  resultado += ")";
      resultado += ";";
   }
   return resultado;
}
	
void TTabConversao::converte(TIndividuo *individuo, TIndividuo *indReduzido)
{
   TGene *gIR;
   TTipoConversao c;
   TGene *gI;
   TGene *gTemp;
   int aux[individuo->get_qtdeGenes()];
   int indG;
   indG = 0;
   for (int i=0; i< indReduzido->get_qtdeGenes(); i++)
   {
      gIR = indReduzido->get_por_indice(i);
	  c = VP_tabConv[gIR->id];
      gI = individuo->get_por_id(c.pr); 

	  do
	  {
         gTemp = gI;
		 aux[indG] = gI->id;
		 gI = gI->prox;
		 indG++;
	  }while(gTemp->id!=c.ul);	  
   }
   
   for (int i = 1; i< individuo->get_qtdeGenes()-1; i++)
   {
      gI = individuo->get_por_indice(i);
      if (gI->id!=aux[i])
	     individuo->troca (gI->id, aux[i]);
   }
}