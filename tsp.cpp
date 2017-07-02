/*
*  tsp.cpp
*
*  Módulo responsável pelo carregamento das instancias a serem testadas
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 05/05/2017
*  Datas de Modificações:
*
*/

#include "tsp.hpp"

//Metodos Privados
int TMapaGenes::getNumGeneDoArquivo(xmlDocPtr doc, xmlNode * a_node)
{
    xmlNode *cur_node = NULL;
    xmlChar *key;
    int val;

    for (cur_node = a_node; cur_node; cur_node = cur_node->next)
    {
       if (cur_node->type == XML_ELEMENT_NODE)
       {
          if (!xmlStrcmp(cur_node->name, (xmlChar *)"description"))
          {
             key = xmlNodeListGetString(doc, cur_node->xmlChildrenNode, 1);
             val = atoi((char *)key);
             xmlFree(key);
             return val;
          }
       }
    }
    return 0;
 }

void TMapaGenes::preencheMapaDist (int geneOri, xmlDocPtr doc, xmlNode * a_node)
{
    xmlNode *cur_node = NULL;
    xmlChar *uri;
    xmlChar *key;
    char dist[25];
    int geneDest;
    float df;

    for (cur_node = a_node; cur_node; cur_node = cur_node->next)
    {
       if (cur_node->type == XML_ELEMENT_NODE)
       {
          uri = xmlGetProp(cur_node, (const xmlChar *)"cost");
          sprintf(dist, "0%s", uri);
          df = atof (dist);

          key = xmlNodeListGetString(doc, cur_node->xmlChildrenNode, 1);
          geneDest = atoi((char *)key);

          set_distancia(geneOri, geneDest, df);

          xmlFree(key);
    xmlFree(uri);
       }

    }
 }

void TMapaGenes::preencheMapa(xmlDocPtr doc, xmlNode * a_node)
{
    xmlNode *cur_node = NULL;
    int gene=0;

    for (cur_node = a_node; cur_node; cur_node = cur_node->next)
    {
       if (cur_node->type == XML_ELEMENT_NODE)
       {
          if (!xmlStrcmp(cur_node->name, (xmlChar *)"vertex"))
          {
             preencheMapaDist(gene++, doc, cur_node->children);
          }
       }
       preencheMapa(doc, cur_node->children);
    }
}

TMapaGenes::TMapaGenes ()
{
      VP_qtdeGenes = -1;
}

TMapaGenes::TMapaGenes (int numGenes)
{
   inicializa (numGenes);
}

int TMapaGenes::get_qtdeGenes () { return VP_qtdeGenes; };

void TMapaGenes::carregaDoArquivo(char *nomeArquivo)
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

       //Obtendo o número de genes no arquivo
       VP_qtdeGenes = getNumGeneDoArquivo(doc, root_element->children);

       //Alocando a tabela
       inicializa (VP_qtdeGenes);

       //preenchendo a tabela com os valores da distáncia
       preencheMapa(doc, root_element->children);

       //liberando documento
       xmlFreeDoc(doc);
       // liberando as variaveis lobais
       xmlCleanupParser();

}

void TMapaGenes::inicializa (int numGenes)
{
       int i;
       int j;
       VP_qtdeGenes = numGenes;
       VP_mapaDist = (double **) malloc(numGenes*sizeof(double *));

       for (i=0; i<VP_qtdeGenes; i++)
       {
          VP_mapaDist [i] = (double *) malloc(numGenes*sizeof(double));
          for (j=0; j<VP_qtdeGenes; j++)
          {
             VP_mapaDist[i][j] = infinito; //Inicia Todos os genes com valor infinito na distância
                                           //ou seja, não tem caminho entre eles
          }
          VP_mapaDist[i][i] = 0.0; //a distância de um gene para ele mesmo é 0
       }

}

TMapaGenes::~TMapaGenes ()
{
       int i;

       for (i=VP_qtdeGenes-1; i>=0; i--)
       {
          free (VP_mapaDist [i]);
       }

       if (VP_qtdeGenes>0) free (VP_mapaDist);

 }

void TMapaGenes::set_distancia(int geneOri, int geneDest, double distancia)
{
       //a distância do gene para ele mesmo não pode ser alterada
       //nenum gene pode está fora do indice d tabela
       if ((geneOri!=geneDest)&&(geneOri>=0)&&(geneOri<VP_qtdeGenes)&&(geneDest>=0)&&(geneDest<VP_qtdeGenes))
          VP_mapaDist[geneOri][geneDest] = distancia;
}

double TMapaGenes::get_distancia(int geneOri, int geneDest)
{
       //nenum gene pode está fora do indice d tabela
       if ((geneOri<VP_qtdeGenes)&&(geneDest>=0)&&(geneDest<VP_qtdeGenes))
          return VP_mapaDist[geneOri][geneDest];
       else
          return 0.0;
}
