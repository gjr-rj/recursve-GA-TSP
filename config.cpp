/*
*  config.cpp
*
*  Módulo responsável pelo carregamento do arquivo de configuração
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 05/05/2017
*  Datas de Modificações:
*
*/

#include "config.hpp"

void TConfig::leInfo(xmlDocPtr doc, xmlNode * a_node)
{
   xmlNode *cur_node = NULL;
   xmlChar *key;
   int val;

   for (cur_node = a_node; cur_node; cur_node = cur_node->next)
   {
      if (cur_node->type == XML_ELEMENT_NODE)
      {
         key = xmlNodeListGetString(doc, cur_node->xmlChildrenNode, 1);
         val = atoi((char *)key);
         xmlFree(key);

         if (!xmlStrcmp(cur_node->name, (xmlChar *)"tamanhoPopulacao")) tamPopulacao = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"numGeracoes")) maxGeracao = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"percentManipulacao")) percentManipulacao = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"percentMutacao")) percentMutacao = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"printParcial")) printParcial = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"percentReducao")) percentReducao = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"numExec")) numExec = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"mutacao")) mutacao = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"cruzamento")) cruzamento = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"percentMutacaoRecursiva")) percentMutacaoRecursiva = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"percentElitismo")) percentElitismo = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"selecao")) selecao = val;
         else if (!xmlStrcmp(cur_node->name, (xmlChar *)"selIndMutacao")) selIndMutacao = val;
		   else if (!xmlStrcmp(cur_node->name, (xmlChar *)"profundidadeMaxima")) profundidadeMaxima = val;

      }

      leInfo(doc, cur_node->children);
   }
}

TConfig::TConfig()
{
   tamPopulacao = 200;
   maxGeracao = 300;
   percentManipulacao = 30;
   percentMutacao = 1;
	percentMutacaoRecursiva = 0;
   percentReducao = 50;
   numExec = 2;
   mutacao = 0;
   cruzamento = 0;
   printParcial = 0;
   percentMutacaoRecursiva = 0;
   percentReducao = 75;
   percentElitismo = 70;
   profundidadeMaxima = 5;
	selecao = 0;
}

void TConfig::carregaDoArquivo(char *nomeArquivo)
{
   xmlDoc *doc = NULL;
   xmlNode *root_element = NULL;

   // Lendo o arquivo
   doc = xmlReadFile(nomeArquivo, NULL, 0);

   if (doc == NULL)
   {
      printf("Erro ao carregar o arquivo %s\n", nomeArquivo);
      return;
   }

   // Obtendo o elemento root
   root_element = xmlDocGetRootElement(doc);

   //preenchendo a tabela com os valores da distáncia
   leInfo(doc, root_element);

   //liberando documento
   xmlFreeDoc(doc);
   // liberando as variaveis lobais
   xmlCleanupParser();
}


