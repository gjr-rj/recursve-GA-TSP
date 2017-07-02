/*
*  tsp.hpp
*
*  Módulo responsável pelo carregamento das instancias a serem testadas
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 05/05/2017
*  Datas de Modificações:
*
*/

#ifndef _TSP_H
#define	_TSP_H


#include "tipos.cpp"
#include <libxml/parser.h>
#include <libxml/tree.h>

/*******************************************************
classe de TMapaGenes. Todas as distâncias entre os genes
********************************************************/
class TMapaGenes
{
   private:
      double **VP_mapaDist;
      int VP_qtdeGenes;

   //Metodos Privados
   int getNumGeneDoArquivo(xmlDocPtr doc, xmlNode * a_node);
   void preencheMapaDist (int geneOri, xmlDocPtr doc, xmlNode * a_node);

   void preencheMapa(xmlDocPtr doc, xmlNode * a_node);

   public:
      TMapaGenes ();

      TMapaGenes (int numGenes);

      int get_qtdeGenes ();

      void carregaDoArquivo(char *nomeArquivo);

      void inicializa (int numGenes);

      ~TMapaGenes ();

      void set_distancia(int geneOri, int geneDest, double distancia);

      double get_distancia(int geneOri, int geneDest);

};

#endif
