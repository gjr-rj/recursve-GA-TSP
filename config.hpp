/*
*  config.hpp
*
*  Módulo responsável pelo carregamento do arquivo de configuração
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 05/05/2017
*  Datas de Modificações:
*
*/

#ifndef _CONFIG_H
#define	_CONFIG_H


#include <libxml/parser.h>
#include <libxml/tree.h>

/*********************************************************
Classe de configuração, utilizada para acelerar o processo
de testes, estava ruim compilar para cada configuração e
estava ruim digitar todos os parâmetros de confguração
**********************************************************/
class TConfig
{
   private:

   void leInfo(xmlDocPtr doc, xmlNode * a_node);

   public:
      /*******************************************
           Parâmetros de teste do sistema
      ********************************************/
      int numExec;              //número de execuções
      int tamPopulacao;
      int maxGeracao;
      int mutacao;                 //mutação base a ser realzada
      int cruzamento;              //cruzamento base a ser realizado
      int percentManipulacao;      //percentual de manipulação do indivíduo. (exclusão / cruzamento)
      int percentMutacao;          //percentual de mutação
	   int profundidadeMaxima;      //Limite para profundidade recursiva
      int printParcial;            //Se inprime informações intermediárias
      int percentMutacaoRecursiva; //percentual de recursividade
      int percentReducao;          //percentual de redução do gene a cada recursvidade
      int percentElitismo;         //Percentual de elitismo
      int selecao;                 //seleção base a ser usada
      int selIndMutacao;           //Seleção do individuos para mutacçao
                                   //  0 - o melhor individuo nunca é selecionado
                                   //  1 - o melhor individuo sempre é seleconado
                                   //  2 - seleção do melhor individuo é aleatória

      /*******************************************************
           Os valores padrões são os utilizados no artigo
      ********************************************************/
      TConfig();
      void carregaDoArquivo(char *nomeArquivo);
};

#endif
