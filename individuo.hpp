/*
*  individuo.hpp
*
*  Módulo responsável pelo controle do individuo
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 13/05/2017
*  Datas de Modificações:
*
*/

#ifndef _INDIVIDUO_H
#define	_INDIVIDUO_H


#include <vector>
#include "tsp.hpp"
#include "arqlog.hpp"
#include "utils.hpp"
#include <algorithm>

/**************************************************
class indivíduo. Uma sequência de genes encadeados
***************************************************/

class TIndividuo
{
   private:
      int VP_Extra;

      vector <TGene *> VP_direto;
      vector <TGene *> VP_indice;

      int VP_qtdeGenes;
      double VP_dist;
		double VP_distInv; // 1/VP_dist para roleta

      int VP_qtdeGeneAlloc;

      TMapaGenes *VP_Mapa;
      TArqLog *VP_ArqSaida;

   public:

      //Apenas leitura
      TMapaGenes *getMapa();
      TArqLog *getArqLog();

      //Propriedades
      int get_extra ();
      void set_extra (int val);

      //Métodos
      TIndividuo (TMapaGenes *mapa, TArqLog *arqSaida);

      ~TIndividuo();

      TGene *get_por_indice (int indice);
      TGene *get_por_id (int id);
      TGene *get_ini ();

      int ind_prox(int ind);
      TGene *prox (TGene *gene);
      TGene *prox (int idG);

      int ind_ant(int ind);
      TGene *ant (TGene *gene);
      TGene *ant (int idG);

      int get_qtdeGenes();

      double get_distancia ();
		double get_distanciaInv ();

      //Cria um novo indivíduo
      void novo ();
      //Cria um novo indivíduo a partir de um vetor de genes
		void novo (vector<TGene *> genes);
	  //Cria um novo individuo a partir da tabela de conversão
      void novo (vector<TTipoConversao> genes, int extraPadrao);

      string toString ();
      string toString (int init);
		string toStringExpand ();

      //Troca a posição de 2 genes pelo nome
      void troca(int g1, int g2);
      //Troca a posição de 2 genes pelo indice
      void troca_indice(int index1, int index2);
	  
		//Troca 2 blocos de genes pelo indice
		//Obs.: Se o tamanho de um dos blocos for 0,
		//      será como se movesse o outro.
		//      Se o tamanho for negativo,
		//      além da troca, haverá uma inversão do bloco        
		void troca_sub(int index1, int size1, int index2, int size2);
	  
      //Embaralha os genes de um individuo
      void embaralha ();
      void embaralha (int index1, int index2);
      //Recalcula o valor de um indivduo
      void recalcDist ();

      //Calcula a distância total no caso de uma sub-rota ser invertida
      double get_dist_sub_reverso(TGene *G1, TGene *G2);
      double get_dist_sub_reverso_id(int id1, int id2);
      double get_dist_sub_reverso_indice(int indice1, int indice2);

      //Inverte o caminho de uma sub-rota
      void inverte_sub_indice (int indice1, int indice2);
      void inverte_sub (TGene *G1, TGene *G2);
      void inverte_sub_id (int id1, int id2);

      //Clona o individuo. Necessário para o elitismo
      TIndividuo *clona ();
      int melhorPossivel ();
};

#endif
