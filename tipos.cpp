/*
*  tipos.hpp
*
*  Módulo responsável pelos típos, constantes e estruturas
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 06/05/2017
*  Datas de Modificações:
*
*/

#ifndef _TIPOS_C
#define	_TIPOS_H

const double infinito = 1000000.00;

/*******************************************************
             Estruturas Básicas
********************************************************/
/*************************************************************
    Estrutura referente a um gene, representando 1 cidade
Contém 2 ponteriros para o anterior e posterior. Também
duas informações que representam a cidade referente quando
está é destino e quando é origem
*************************************************************/
struct TGene
{
   int id;
   int ori;   //O gene sendo a origem
   int dest;  //O gene sendo o detstino
   int i;     //A posição de um gene
   struct TGene *prox;
   struct TGene *ant;
};

struct TTipoConversao
{
   int pr;
   int prOrig;
   int prDest;
   int ul;
   int ulOrig;
   int ulDest;
};

#endif
