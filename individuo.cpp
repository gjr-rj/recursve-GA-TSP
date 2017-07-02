/*
*  individuo.cpp
*
*  Módulo responsável pelo controle do individuo
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 13/05/2017
*  Datas de Modificações:
*
*/

#include "individuo.hpp"
#include <iostream>

//Apenas letura
TMapaGenes *TIndividuo::getMapa() { return VP_Mapa; }
TArqLog *TIndividuo::getArqLog()  { return VP_ArqSaida; }

int TIndividuo::get_extra ()         { return VP_Extra; }
void TIndividuo::set_extra (int val) { VP_Extra = val; }

double TIndividuo::get_distancia() { return VP_dist; }

TGene *TIndividuo::get_ini () { return VP_indice[0]; }
TGene *TIndividuo::get_por_indice (int indice) { return VP_indice[indice]; }
TGene *TIndividuo::get_por_id (int id) { return VP_direto[id]; }

int TIndividuo::ind_prox(int ind) { return (ind<(VP_qtdeGenes-1))?ind+1:0; }
TGene *TIndividuo::prox (TGene *gene) { return gene->prox; }
TGene *TIndividuo::prox (int idG) { return VP_direto[idG]->prox; }

int TIndividuo::ind_ant(int ind) { return (ind>0)?ind-1:(VP_qtdeGenes-1); }
TGene *TIndividuo::ant (TGene *gene) { return gene->ant; }
TGene *TIndividuo::ant (int idG) { return VP_direto[idG]->ant; }

int TIndividuo::get_qtdeGenes() { return VP_qtdeGenes; }

//Métodos
TIndividuo::TIndividuo (TMapaGenes *mapa, TArqLog *arqSaida)
{
   VP_Mapa = mapa;
   VP_ArqSaida = arqSaida;

   VP_Extra = 0;
   VP_dist = 0;
   VP_qtdeGeneAlloc = 0;
}

TIndividuo::~TIndividuo()
{
   for (vector< TGene *>::iterator i=VP_indice.begin(); i!=VP_indice.end(); ++i)
   {
     delete (*i);
   }
   VP_indice.clear();
   VP_direto.clear();
}

//Cria um novo indivíduo
void TIndividuo::novo ()
{
   TGene *gene;
   VP_qtdeGenes = VP_Mapa->get_qtdeGenes();
   VP_dist = 0;

   for (int i=0; i<VP_qtdeGenes; i++)
   {
      gene = new TGene;
      gene->id   = i;
      gene->ori  = i;
      gene->dest = i;
      gene->i    = i;
      VP_direto.push_back (gene);
      VP_indice.push_back (gene);
      if (i==0) continue;
      VP_direto[i]->ant = VP_direto[i-1];
      VP_direto[i-1]->prox = VP_direto[i];
      VP_dist += VP_Mapa->get_distancia(VP_direto[i-1]->ori, VP_direto[i]->dest);
   }
   //Fechando o ciclo
   VP_direto[VP_qtdeGenes-1]->prox = VP_direto[0];
   VP_direto[0]->ant = VP_direto[VP_qtdeGenes-1];

   //Se o indivíduo foi modificado, o 2opt deve começar do início
   VP_Extra = 0;

   VP_dist += VP_Mapa->get_distancia(VP_direto[VP_qtdeGenes-1]->ori, VP_direto[0]->dest);
}

void TIndividuo::novo (vector<TTipoConversao> genes, int extraPadrao)
{
   TGene *gene;
   VP_qtdeGenes = genes.size();
   VP_dist = 0;

   for (int i=0; i<VP_qtdeGenes; i++)
   {
      gene = new TGene;
      gene->id   = i;
      gene->ori  = genes[i].ulOrig;
      gene->dest = genes[i].prDest;
      gene->i    = i;

      VP_direto.push_back (gene);
      VP_indice.push_back (gene);

      if (i==0) continue;
      VP_indice[i]->ant = VP_indice[i-1];
      VP_indice[i-1]->prox = VP_indice[i];

      VP_dist += VP_Mapa->get_distancia(VP_indice[i-1]->ori, VP_indice[i]->dest);
   }

   //Fechando o ciclo
   VP_indice[VP_qtdeGenes-1]->prox = VP_indice[0];
   VP_indice[0]->ant = VP_indice[VP_qtdeGenes-1];

   //Se o indivíduo foi modificado, o 2opt deve começar do início
   VP_Extra = extraPadrao;

   VP_dist += VP_Mapa->get_distancia(VP_indice[VP_qtdeGenes-1]->ori, VP_indice[0]->dest);	
}

//Cria um novo indivíduo de um vetor de TGene
void TIndividuo::novo (vector<TGene *> genes)
{
   TGene *gene;
   VP_qtdeGenes = genes.size();
   VP_direto.resize (VP_qtdeGenes);
   VP_dist = 0;

   for (int i=0; i<VP_qtdeGenes; i++)
   {
      gene = new TGene;
      gene->id   = genes[i]->id;
      gene->ori  = genes[i]->ori;
      gene->dest = genes[i]->dest;
      gene->i    = i;

      VP_direto[gene->id] = gene;
      VP_indice.push_back (gene);

      if (i==0) continue;
      VP_indice[i]->ant = VP_indice[i-1];
      VP_indice[i-1]->prox = VP_indice[i];

      VP_dist += VP_Mapa->get_distancia(VP_indice[i-1]->ori, VP_indice[i]->dest);
   }

   //Fechando o ciclo
   VP_indice[VP_qtdeGenes-1]->prox = VP_indice[0];
   VP_indice[0]->ant = VP_indice[VP_qtdeGenes-1];

   //Se o indivíduo foi modificado, o 2opt deve começar do início
   VP_Extra = 0;

   VP_dist += VP_Mapa->get_distancia(VP_indice[VP_qtdeGenes-1]->ori, VP_indice[0]->dest);
}

string TIndividuo::toString ()
{
   string resultado = "";
   for (vector< TGene *>::iterator i=VP_indice.begin(); i!=VP_indice.end(); ++i)
   {
     resultado += to_string( (*i)->id);
     resultado += ";";
   }
   return resultado;
}

string TIndividuo::toStringExpand ()
{
   string resultado = "";
   for (vector< TGene *>::iterator i=VP_indice.begin(); i!=VP_indice.end(); ++i)
   {
     resultado += to_string( (*i)->id);
	 resultado += "[";
     resultado += to_string( (*i)->ori );
     resultado += ",";
     resultado += to_string( (*i)->dest );
     resultado += "] ";
     resultado += "(";
     resultado += to_string( (*i)->ant->id );
     resultado += ",";
     resultado += to_string( (*i)->prox->id );
     resultado += " / ";
     resultado += to_string( (*i)->i);
     resultado += ")";
     resultado += ";";
   }
   return resultado;
}

string TIndividuo::toString (int init)
{
   string resultado = "";
   for (vector< TGene *>::iterator i=VP_indice.begin(); i!=VP_indice.end(); ++i)
   {
     resultado += to_string( (*i)->id+1);
     resultado += ";";
   }
   return resultado;
}

void TIndividuo::troca(int g1, int g2)
{
   troca_indice(VP_direto[g1]->i, VP_direto[g2]->i);
}

void TIndividuo::troca_indice(int index1, int index2)
{
   if ((index1==0)||(index2==0)||(index1==index2)) return;

   TGene *tempG;

   if(index1>index2)
   {
      int i = index1;
      index1 = index2;
      index2 = i;
   }

   VP_dist -= VP_Mapa->get_distancia(VP_indice[ind_ant(index1)]->ori, VP_indice[index1]->dest);
   VP_dist -= VP_Mapa->get_distancia(VP_indice[index2]->ori, VP_indice[ind_prox(index2)]->dest);
   VP_dist -= VP_Mapa->get_distancia(VP_indice[index1]->ori, VP_indice[ind_prox(index1)]->dest);

   //se = distância de 1 a 1+1 = distância de 2-1 a 2
   //Seria contado duas vezes
   if((ind_prox(index1))!=index2)
      VP_dist -= VP_Mapa->get_distancia(VP_indice[ind_ant(index2)]->ori, VP_indice[index2]->dest);


   tempG = VP_indice[index1];
   VP_indice[index1] = VP_indice[index2];
   VP_indice[index2] = tempG;
   VP_indice[index1]->i = index1;
   VP_indice[index2]->i = index2;

   //resolvendo o encadeamento para index1
   VP_indice[ind_ant(index1)]->prox = VP_indice[index1];
   VP_indice[index1]->ant = VP_indice[ind_ant(index1)];
   VP_indice[index1]->prox = VP_indice[ind_prox(index1)];
   VP_indice[ind_prox(index1)]->ant = VP_indice[index1];

   //resolvendo o encadeamento para index2
   VP_indice[ind_ant(index2)]->prox = VP_indice[index2];
   VP_indice[index2]->ant = VP_indice[ind_ant(index2)];
   VP_indice[index2]->prox = VP_indice[ind_prox(index2)];
   VP_indice[ind_prox(index2)]->ant = VP_indice[index2];

   //Somando as novas dstâncias
   VP_dist += VP_Mapa->get_distancia(VP_indice[ind_ant(index1)]->ori, VP_indice[index1]->dest);
   VP_dist += VP_Mapa->get_distancia(VP_indice[index2]->ori, VP_indice[ind_prox(index2)]->dest);
   VP_dist += VP_Mapa->get_distancia(VP_indice[index1]->ori, VP_indice[ind_prox(index1)]->dest);

   //se = distância de 1 a 1+1 = distância de 2-1 a 2
   //Seria contado duas vezes
   if((ind_prox(index1))!=index2)
      VP_dist += VP_Mapa->get_distancia(VP_indice[ind_ant(index2)]->ori, VP_indice[index2]->dest);

}

/*
	  Troca 2 blocos de genes, pelo índice.
	  Obs.: Se o tamanho de um dos blocos for 0,
	        será como se movesse o outro.
	        Se o tamanho for negativo,
	        além da troca, haverá uma inversão do bloco        
*/
void TIndividuo::troca_sub(int index1, int size1, int index2, int size2)
{
   int ini1, ini2, increment1, increment2;

   //Por convenção, a rota sempre iniciará no gene 0
   //logo, este não pode ser trocado
   if (index1==0) return;
   
   //Não há troca a ser feita
   if ((size1==size2)&&(size1==0)) return;

   if (size1>0) 
   { 
      ini1 = index1;
      increment1 = 1;
   }
   else //Inverte o bloco
   {
      size1 *= -1;
      ini1 = index1+size1-1;
      increment1 = -1;
   }

   if (size2>0) 
   { 
      ini2 = index2;
      increment2 = 1;
   }
   else //Inverte inverte o bloco
   {
      size2 *= -1;
      ini2 = index2+size2-1;
      increment2 = -1;
   }
   
   //Não pode sobrepor
   if (index1+size1>index2) return;
   //Não pode ultrapassar
   if (index2+size2>this->get_qtdeGenes()) return;
 
   //Montando os bolocos
   int v1[size1];
   int v2[size2];
   int vm[index2-index1-size1];
  
   //Montando o prmeiro bloco
   for (int i = 0, j=ini1; i<size1; i++, j+=increment1) v1[i] = this->VP_indice[j]->id; 
   //Montando a parte entre os 2 blocos
   for(int i = 0, j=index1+size1; i<index2-index1-size1; i++, j++) vm[i] = this->VP_indice[j]->id;
   //Montando o segundo bloco
   for (int i = 0, j=ini2; i<size2; i++, j+=increment2) v2[i] = this->VP_indice[j]->id;

   //Remontando o indivíduo
   //o segundo bloco passa a ser o prmeiro
   for (int i = 0, j=index1; i<size2; i++, j++) this->troca (this->VP_indice[j]->id, v2[i]); 
   //O bloco do meio
   for(int i = 0, j=index1+size2; i<index2-index1-size1; i++, j++) this->troca (this->VP_indice[j]->id,vm[i]);
   //o prmeiro bloco passará a ser o segundo
   for (int i = 0, j=index2+size2-size1; i<size1; i++, j++) this->troca (this->VP_indice[j]->id, v1[i]); 	
}

void TIndividuo::embaralha ()
{
   embaralha(1, get_qtdeGenes()-1);
}

void TIndividuo::embaralha (int index1, int index2)
{
   int rd1;
   int rd2;

   if (index1>index2)
   {
      rd1 = index1;
      index1 = index2;
      index2 = rd1;
   }

   for (int i=0; i<=((index2-index1)/2); i++)
   {
      rd1 = TUtils::rnd(index1, index2);
      rd2 = TUtils::rnd(index1, index2);
      troca_indice(rd1, rd2);
   }
}

void TIndividuo::recalcDist ()
{
   VP_dist = 0;

    for (int i = 0; i<VP_qtdeGenes; i++)
       VP_dist += VP_Mapa->get_distancia(VP_indice[ind_ant(i)]->ori, VP_indice[i]->dest);
}

double TIndividuo::get_dist_sub_reverso(TGene *G1, TGene *G2)
{
   TGene *tempG;

   if ((G1->i==0)||(G2->i==0)||(G1->i == G2->i)) return VP_dist;
   else if (G1->i > G2->i)
   {
      tempG = G1;
      G1 = G2;
      G2 = tempG;
   }

   double tot=VP_dist;

   //Entre os genes 1 e 2, a soma é do caminho de retorno ou seja, de 2 para 1
   for (tempG=G1; tempG->id!=G2->id; tempG = tempG->prox)
   {
      tot -= VP_Mapa->get_distancia(tempG->ori, tempG->prox->dest);
      tot += VP_Mapa->get_distancia(tempG->prox->ori, tempG->dest);
   }

   //Arestas que ficaram faltando
   tot -= VP_Mapa->get_distancia(G1->ant->ori, G1->dest);
   tot -= VP_Mapa->get_distancia(G2->ori, G2->prox->dest);

   tot += VP_Mapa->get_distancia(G1->ant->ori, G2->dest);
   tot += VP_Mapa->get_distancia(G1->ori, G2->prox->dest);

   return tot;
}

double TIndividuo::get_dist_sub_reverso_id(int id1, int id2)
{
   return get_dist_sub_reverso (VP_direto[id1], VP_direto[id2]);
}

double TIndividuo::get_dist_sub_reverso_indice(int indice1, int indice2)
{
   return get_dist_sub_reverso (VP_indice[indice1], VP_indice[indice2]);
}

void TIndividuo::inverte_sub_indice (int indice1, int indice2)
{
   int i;
   if ((indice1==0)||(indice2==0)||(indice1==indice2)) return;
   else if (indice1>indice2)
   {
      i = indice1;
      indice1 = indice2;
      indice2 = i;
   }

   for (i=0; i<=(indice2-indice1)/2;i++)
      troca_indice(indice1+i, indice2-i);
}

void TIndividuo::inverte_sub (TGene *G1, TGene *G2)
{
   inverte_sub_indice (G1->i, G2->i);
}

void TIndividuo::inverte_sub_id (int id1, int id2)
{
   inverte_sub_indice (VP_direto[id1]->i, VP_direto[id2]->i);
}

//Clona o individuo. Necessário para o elitismo
TIndividuo *TIndividuo::clona ()
{
   TIndividuo *clone = new TIndividuo(getMapa(), getArqLog());
   clone->novo(VP_indice);
   clone->set_extra(get_extra());
   return clone;
}

int TIndividuo::melhorPossivel()
{
   int melhorou = 0; //retrna false (=0) quando não mehora ou true (!=0) quando melhora 
   vector<TGene *> melhor;
   vector<TGene *> temp;
   double melhorDist = get_distancia();
   double distTemp;

   //Realizando todas as combinações, fixando o 0, pois é um ciclo
   for(int i=1;i<get_qtdeGenes();++i) temp.push_back(get_por_id(i));

   do
   {
      //Calculando o peso total de cada permutação
      distTemp = VP_Mapa->get_distancia(VP_direto[0]->ori, temp[0]->dest);
      for (unsigned i=1;i<temp.size();++i) distTemp += VP_Mapa->get_distancia(temp[i-1]->ori, temp[i]->dest);
      distTemp += VP_Mapa->get_distancia(temp[temp.size()-1]->ori ,VP_direto[0]->dest);

      if (distTemp<melhorDist)
      {
         melhorDist = distTemp;
         melhor = temp;
      }
   }while(next_permutation(temp.begin(),temp.end()));

   //Se encontrei combinação melhor do que a existente
   if (melhorDist<get_distancia())
   {
      melhorou = 1;
      //Mudo o individuo para mesma sequência do melhor
      for (int i=1;i<get_qtdeGenes()-1;++i)
      {
         if (VP_indice[i]->id!=melhor[i-1]->id)
            troca_indice(i, (get_por_id(melhor[i-1]->id))->i);
      }
   }
   
   melhor.clear();
   temp.clear();
   
   return melhorou;
}

