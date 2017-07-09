/*
*  populacao.cpp
*
*  Módulo responsável tratamento da população
*  PPGI - NCE - UFRJ
*  Data Criação: 13/05/2017
*  Datas de Modificações:
*
*/

#include "populacao.hpp"
#include <iostream>

//Apenas leitura
TMapaGenes *TPopulacao::getMapa()           { return VP_Mapa; }
TArqLog *TPopulacao::getArqLog()            { return VP_ArqSaida; }
unsigned TPopulacao::get_tamanho()          { return VP_tamanho; }
unsigned TPopulacao::get_qtdeIndividuo()    { return VP_individuos.size(); }
double TPopulacao::get_soma_dist ()         { return VP_somaDistancias; }
double TPopulacao::get_soma_distInv ()      { return VP_somaDistanciasInv; }

//Métodos
TIndividuo *TPopulacao::get_individuo(unsigned indice)
{
   if (indice<0) indice = 0;
   if (indice>=VP_tamanho) indice = VP_tamanho-1;

   return VP_individuos[indice];
}

void TPopulacao::sub_dist_tot (double val)  
{ 
	VP_somaDistancias -= val; 
	VP_somaDistanciasInv -= 1/val;
}
void TPopulacao::soma_dist_tot (double val) 
{ 
	VP_somaDistancias += val; 
	VP_somaDistanciasInv += 1/val;
}

void TPopulacao::add_individuo(TIndividuo *individuo)
{
   if (VP_individuos.size() < VP_tamanho)
   {
      VP_individuos.push_back(individuo);
      VP_somaDistancias += individuo->get_distancia();
		VP_somaDistanciasInv += individuo->get_distanciaInv();
   }
   else
	   delete individuo;
}

void TPopulacao::add_individuo(vector <TIndividuo *>individuos)
{
   for (unsigned i =0; i < individuos.size(); i++ )
      add_individuo(individuos[i]);
		
  individuos.clear(); 
}

TPopulacao::TPopulacao (int tamanho, TMapaGenes *mapa, TArqLog *arqSaida)
{
   VP_Mapa = mapa;
   VP_ArqSaida = arqSaida;

   VP_somaDistancias = 0;
	VP_somaDistanciasInv = 0;
   VP_tamanho = tamanho;
}

void TPopulacao::povoa()
{
   TIndividuo *individuo;

   for (unsigned i=0; i<VP_tamanho; i++)
   {
      individuo = new TIndividuo(getMapa(), getArqLog());
      individuo->novo();
      individuo->embaralha();
      VP_somaDistancias += individuo->get_distancia();
		VP_somaDistanciasInv += individuo->get_distanciaInv();
      VP_individuos.push_back(individuo);
   }
}

void TPopulacao::povoa(vector<TTipoConversao> tabConv, int extraPadrao)
{
   TIndividuo *individuo;

   for (unsigned i=0; i<VP_tamanho; i++)
   {
      individuo = new TIndividuo(getMapa(), getArqLog());
      individuo->novo(tabConv, extraPadrao);
	  
	  //Não embaralha o primeiro individuo
	  //Pois a sequencia 1,2,3, ..., n, n+1
	  //é o melhor individuo do nivel anterior
	  if (i)  individuo->embaralha();
	  
      VP_somaDistancias += individuo->get_distancia();
		VP_somaDistanciasInv += individuo->get_distanciaInv();
      VP_individuos.push_back(individuo);
   }
}

TPopulacao::~TPopulacao()
{
   for (vector< TIndividuo *>::iterator i=VP_individuos.begin(); i!=VP_individuos.end(); ++i)
   {
     delete (*i);
   }
   VP_individuos.clear();
}

string TPopulacao::toString ()
{
   string resultado = "";
   for (vector< TIndividuo *>::iterator i=VP_individuos.begin(); i!=VP_individuos.end(); ++i)
   {
     resultado += (*i)->toString();
     resultado += ";(";
     resultado += to_string((*i)->get_distancia());
     resultado += ")";
     resultado += "\n";
   }
   return resultado;
}

void TPopulacao::ordena()
{
   //sort (VP_individuos.begin(), VP_individuos.end(), comparaMelhor);
   qsort (&VP_individuos[0], VP_individuos.size(), sizeof (TIndividuo *), comparaMelhor);
}

int TPopulacao::comparaMelhor(const void *Indi, const void *Indj)
{
   //As linhas comentadas referen-se a comparação com o sort
   //TIndividuo *i =  (TIndividuo *) Indi;
   //TIndividuo *j =  (TIndividuo *) Indj;

   TIndividuo *i =  *(TIndividuo **) Indi;
   TIndividuo *j =  *(TIndividuo **) Indj;

   //return i->get_distancia() < j->get_distancia();
   if(i->get_distancia() > j->get_distancia()) return 1;
   if(i->get_distancia() < j->get_distancia()) return -1;
   return 0;
}

void TPopulacao::troca (int indice1, int indice2)
{
   TIndividuo *temp;
   temp = VP_individuos[indice1];
   VP_individuos[indice1] = VP_individuos[indice2];
   VP_individuos[indice2] = temp;
}

TIndividuo *TPopulacao::get_melhor() { return VP_individuos[0]; }
TIndividuo *TPopulacao::get_pior()   { return VP_individuos[VP_tamanho-1]; }
double TPopulacao::distanciaMedia () { return VP_somaDistancias/VP_tamanho; }
