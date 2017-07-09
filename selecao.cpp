/*
*  selecao.cpp
*
*  Módulo responsável pelas reotinas de seleção
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 08/07/2017
*  Datas de Modificações:
*
*/

#include "selecao.hpp"
#include <iostream>

//Apenas letura
TMapaGenes *TSelecao::getMapa() { return VP_Mapa; }
TArqLog *TSelecao::getArqLog()  { return VP_ArqSaida; }

int TSelecao::getTipo()  { return VP_TipoSelecao; }
void TSelecao::setTipo(int val)  { VP_TipoSelecao = val; }

TSelecao::TSelecao (TMapaGenes *mapa, TArqLog *arqSaida, int tipoSelecao)
{
   VP_Mapa = mapa;
   VP_ArqSaida = arqSaida;
	VP_TipoSelecao = tipoSelecao;
}

TSelecao::~TSelecao () { }

unsigned TSelecao::processa(TPopulacao *populacao, unsigned exceto, int k)
{
	unsigned selecionado = 0;
	switch (this->VP_TipoSelecao)
   {
      case 0:
      {
         selecionado = sorteio (populacao);
         break;
      }
      case 1:
      {
         selecionado = torneioK (populacao, k);
         break;
      }
      case 2:
      {
			selecionado = roleta (populacao);
         break;
      }
      case 3:
      {
			selecionado = roletaInv (populacao);
         break;
      }
      default:
      {
         break;
      }
   }

   if (selecionado==exceto) 
	{
		if (selecionado==0) selecionado++;
		else                selecionado--;
   }

   return selecionado;	
}

unsigned TSelecao::processa(TPopulacao *populacao, unsigned exceto)
{
   return this->processa(populacao, exceto, 5);	
}

unsigned TSelecao::processa(TPopulacao *populacao)
{
   return this->processa(populacao, -1, 5);
}

/**
 *
 * Selección de Torneo K
 *
 **/
unsigned TSelecao::torneioK(TPopulacao *populacao, int k)
{
	unsigned sel = populacao->get_tamanho() - 1;
	TIndividuo *ind = populacao->get_individuo(sel);
   double valSel = ind->get_distancia();
	int opc;

	for(int i = 0; i < k; i++)
	{
		opc = TUtils::rnd(0, populacao->get_tamanho() - 1);
		ind = populacao->get_individuo(opc);
		if(valSel > ind->get_distancia())
		{
			sel = opc;
			ind = populacao->get_individuo(sel);
			valSel = ind->get_distancia();			
		}
  }

  return sel;
}

unsigned TSelecao::sorteio(TPopulacao *populacao)
{
   return TUtils::rnd(0, populacao->get_tamanho() - 1);
}

unsigned TSelecao::roleta(TPopulacao *populacao)
{
	double pick = TUtils::rndd(0, populacao->get_soma_distInv());
	double sum = 0.0;
  
	unsigned i;
	TIndividuo *Ind;
	for(i = 0; (sum < pick) && (i < populacao->get_tamanho()); i++)
	{
		Ind = populacao->get_individuo(i);
      sum += Ind->get_distanciaInv();
	}

	return i-1;
}

unsigned TSelecao::roletaInv(TPopulacao *populacao)
{
	double pick = TUtils::rnd(0, populacao->get_soma_dist());
	double sum = 0.0;
  
	unsigned i;
	TIndividuo *Ind;
	for(i = 0; (sum < pick) && (i < populacao->get_tamanho()); i++)
	{
		Ind = populacao->get_individuo(i);
      sum += Ind->get_distancia();
	}
  
	return i-1;
}
