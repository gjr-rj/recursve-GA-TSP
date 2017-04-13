/**
 * Autores: Geraldo José Ferreira Chagas Júnior
 *
 *
 * Data: 01/2017
 *
 * Programa de Algorítmo Genetico para resolver tsp de forma recursiva
 * Necessário instalação da biblioteca libxml
 * Instalação da biblioteca no debian
 *           $> sudo apt-get install libxml
 *
 * compilação:
 *           $> g++ `xml2-config --cflags --libs` -o tsp tsp.cpp
 *      ou
 *           $> g++ -Wall -pedantic `xml2-config --cflags --libs` -o tsp tsp.cpp
 *     se desejar ver os warnings
 *
 *
 * Testes:
 *      1) Mutação 2opt, apenas se melhorar a rota X Mutção 2opt sempre 
 *                                                   (Sempre multacção no melhor X aleatório para o melhor)
 *      2) Nova geração substituindo os piores     X Nova geração substituindo a menor diferença
 *      3) Sorteio de pares, todos o mesmo peso    X Sorteio de pares com peso
 *      4) cross over de pares aleatórios          X Sorteio de pares p1 e p2 em ordem
 *                                                   (melhor prmeiro X pior primeiro)
 */
#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef LIBXML_TREE_ENABLED

const double infinito = 1000000.00;

int compare(const void *x, const void *y);

time_t sysTime1, sysTime2;

/*******************************************************
             Estruturas Básicas
********************************************************/
/*************************************************************
    Estrutura referente a um gene, representando 1 cidade
Contém 2 ponteriros para o que seria o gene anterior
e posterior.
*************************************************************/
struct TGene
{
   int id;
   int ori;  //O gene sendo a origem
   int dest; //O gene sendo o detstino
   struct TGene *prox;
   struct TGene *ant;   
};

struct TTipoConvercao
{
   int id;
   int qtde;
   int *gene;
   int ori;  //Qual gene original da matriz de distância servirá de orgem para o gene
   int dest; //Qual gene original da matriz de distância servirá de destino para o gene
   struct TTipoConvercao *prox;
   struct TTipoConvercao *ant;
};


/************************************************************************************
No lugar de utilizar uma função de n multiplcações em um algorítmo genético,
será utilizado a aproximação de Stirling para agilizar o algoritmo, pois
neste caso a aproximação será satisfatóra e utilizando pouco esforço computacional.
************************************************************************************/
double fatorialStirling (int n)
{
   const double pi = 3.14;
   const double e  = 2.72;

   long result = 1;
 
   //A aproximação de Stirling funciona apenas para n>6
   if (n==2) result = 2;
   else if (n==3) result = 6;
   else if (n==4) result = 24;
   else if (n==5) result = 120;
   else if (n==6) result = 720;
   else if (n>6)  result = sqrt(2*pi*n)*pow(n/e,n);

   return result;
}

/*******************************************************
classe de TMapaGenes. Todas as distâncias entre os genes
********************************************************/
class TMapaGenes
{
   private:
      double **VP_mapaDist;
      int VP_qtdeGenes;

   //Metodos Privados
   int getNumGeneDoArquivo(xmlDocPtr doc, xmlNode * a_node)
   {
      xmlNode *cur_node = NULL;
      xmlChar *key;
      int val;

      srand(time(NULL));
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

   void preencheMapaDist (int geneOri, xmlDocPtr doc, xmlNode * a_node)
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

   void preencheMapa(xmlDocPtr doc, xmlNode * a_node)
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

   public:
      TMapaGenes ()
      {
        VP_qtdeGenes = -1;
      }

      TMapaGenes (int numGenes)
      {
         inicializa (numGenes);
      }

      int get_qtdeGenes () { return VP_qtdeGenes; };

      void carregaDoArquivo(char *nomeArquivo)
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

      void inicializa (int numGenes) 
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
 
      ~TMapaGenes () 
      {
         int i;

         for (i=VP_qtdeGenes-1; i>=0; i--)
         {
            free (VP_mapaDist [i]);
         }

         if (VP_qtdeGenes>0) free (VP_mapaDist);
      }

      void set_distancia(int geneOri, int geneDest, double distancia)
      { 
         //a distância do gene para ele mesmo não pode ser alterada
         //nenum gene pode está fora do indice d tabela
         if ((geneOri!=geneDest)&&(geneOri>=0)&&(geneOri<VP_qtdeGenes)&&(geneDest>=0)&&(geneDest<VP_qtdeGenes))
            VP_mapaDist[geneOri][geneDest] = distancia;         
      }

      double get_distancia(int geneOri, int geneDest)
      {
         //nenum gene pode está fora do indice d tabela
         if ((geneOri<VP_qtdeGenes)&&(geneDest>=0)&&(geneDest<VP_qtdeGenes))
            return VP_mapaDist[geneOri][geneDest];         
         else
            return 0.0;
      }

};


/************************************************************************
   Fim da clase referente ao carregamento da tabela de distâncias
                Início do Algoritmo Genético
***********************************************************************/

/**************************************************
class indivíduo. Uma sequência de genes encadeados
***************************************************/

class TIndividuo
{
   private:
      int VP_ultOpt;

      TGene **VP_direto;
      TGene *VP_geneIni;
      int VP_qtdeGenes;
      double VP_dist;

      int VP_qtdeGeneAlloc;

   /**************************************
    funções privadas (métodos privados)
   ***************************************/
   double calcDistTotTroca (int idG1, int idG2)
   {
      TGene *tempG;

      double tot=VP_dist;

      //Entre os genes 1 e 2, a soma é do caminho de retorno ou seja, de 2 para 1
      for (tempG=VP_direto[idG1]; tempG->id!=idG2; tempG = tempG->prox)
      {
         tot -= Mapa->get_distancia(tempG->ori, tempG->prox->dest);
         tot += Mapa->get_distancia(tempG->prox->ori, tempG->dest); 
      }

      //Arestas que ficaram faltando
      tot -= Mapa->get_distancia(VP_direto[idG1]->ant->ori, VP_direto[idG1]->dest);
      tot -= Mapa->get_distancia(VP_direto[idG2]->ori, VP_direto[idG2]->prox->dest);

      tot += Mapa->get_distancia(VP_direto[idG1]->ant->ori, VP_direto[idG2]->dest);
      tot += Mapa->get_distancia(VP_direto[idG1]->ori, VP_direto[idG2]->prox->dest);

      return tot;
   }

   //Troca a posição de dois genes 
   void trocaSimples(int g1, int g2)
   {
      TGene *tempG;

      if(VP_direto[g1]->prox->id==VP_direto[g2]->id)
      {
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->ant->ori, VP_direto[g1]->dest);
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->ori, VP_direto[g1]->prox->dest);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->ori, VP_direto[g2]->prox->dest);

         VP_dist += Mapa->get_distancia(VP_direto[g1]->ant->ori, VP_direto[g2]->dest);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->ori, VP_direto[g1]->dest);
         VP_dist += Mapa->get_distancia(VP_direto[g1]->ori, VP_direto[g2]->prox->dest);

         VP_direto[g1]->ant->prox = VP_direto[g2];
         VP_direto[g2]->ant = VP_direto[g1]->ant;
         VP_direto[g1]->prox = VP_direto[g2]->prox;
         VP_direto[g1]->prox->ant = VP_direto[g1];
         VP_direto[g1]->ant = VP_direto[g2];
         VP_direto[g2]->prox = VP_direto[g1];
      }
      else if (VP_direto[g2]->prox->id==VP_direto[g1]->id)
      {
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->ant->ori, VP_direto[g2]->dest);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->ori, VP_direto[g2]->prox->dest);
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->ori, VP_direto[g1]->prox->dest);

         VP_dist += Mapa->get_distancia(VP_direto[g2]->ant->ori, VP_direto[g1]->dest);
         VP_dist += Mapa->get_distancia(VP_direto[g1]->ori, VP_direto[g2]->dest);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->ori, VP_direto[g1]->prox->dest);

         VP_direto[g2]->ant->prox = VP_direto[g1];
         VP_direto[g1]->ant = VP_direto[g2]->ant;
         VP_direto[g2]->prox = VP_direto[g1]->prox;
         VP_direto[g2]->prox->ant = VP_direto[g2];
         VP_direto[g2]->ant = VP_direto[g1];
         VP_direto[g1]->prox = VP_direto[g2];
      }
      else if (g1!=g2)
      {
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->ant->ori, VP_direto[g1]->dest);
         VP_dist -= Mapa->get_distancia(VP_direto[g1]->ori, VP_direto[g1]->prox->dest);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->ant->ori, VP_direto[g2]->dest);
         VP_dist -= Mapa->get_distancia(VP_direto[g2]->ori, VP_direto[g2]->prox->dest);

         VP_dist += Mapa->get_distancia(VP_direto[g1]->ant->ori, VP_direto[g2]->dest);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->ori, VP_direto[g1]->prox->dest);
         VP_dist += Mapa->get_distancia(VP_direto[g2]->ant->ori, VP_direto[g1]->dest);
         VP_dist += Mapa->get_distancia(VP_direto[g1]->ori, VP_direto[g2]->prox->dest);

         VP_direto[g1]->ant->prox = VP_direto[g2];
         VP_direto[g1]->prox->ant = VP_direto[g2];

         VP_direto[g2]->ant->prox = VP_direto[g1];
         VP_direto[g2]->prox->ant = VP_direto[g1];

         tempG = VP_direto[g1]->prox;
         VP_direto[g1]->prox = VP_direto[g2]->prox;
         VP_direto[g2]->prox = tempG;

         tempG = VP_direto[g1]->ant;
         VP_direto[g1]->ant = VP_direto[g2]->ant;
         VP_direto[g2]->ant = tempG;
      }

      //checa((char *)"trocaSimples");
   }

   // Combina de form recursiva todas as possibilidades de forma a encontrar
   // o melhor indivíduo possível. (fotial de quantidade de genes)
   void combina (TGene *ga, int *melhorDeTodos, double *distMelhor)
   {
      TGene *tempG;
      int gaI;
      int TempGI;
      
      if (ga->prox->id==VP_geneIni->id)
      {
         if (VP_dist<*distMelhor)
         {
	    *distMelhor = VP_dist;
            toArray (melhorDeTodos);
            VP_ultOpt = VP_geneIni->prox->id;
         }
         return;
      }

      for (tempG = ga; tempG->id!=VP_geneIni->id; tempG = tempG->prox)
      {
        gaI = ga->id;
	TempGI = tempG->id;

        //Realiza a troca do gene
	trocaSimples (gaI, TempGI);

        combina (VP_direto[TempGI]->prox, melhorDeTodos, distMelhor);

        //Destroca os genes
        trocaSimples (gaI, TempGI);

        //Como é trabahado com ponteiros, as trocas fazem com que
        //os apontamentos saiam da ordem desejada. A utilização
        //do id antes da troca serve para recalibrar os ponteiros
	tempG = VP_direto[TempGI];
	ga = VP_direto[gaI];
      }
      
      //checa((char *)"combina");
   }

   //O 2opt
   void troca (int g1, int g2)
   {
      if ((g1==0)||(g2==0)||(g1==g2)) return; //A primeira posição nunca será trocada.

      TGene *tempGene; 
      int vizinho;
      
      bool troca = false;
      int guarda;

      VP_dist = 0;

      if (VP_direto[g1]->prox->id==g2)
      {
         int AG1 = VP_direto[g1]->ant->id;
         int PG2 = VP_direto[g2]->prox->id;

          VP_direto[AG1]->prox = VP_direto[g2]; VP_direto[g2]->ant = VP_direto[AG1];
          VP_direto[g2]->prox  = VP_direto[g1]; VP_direto[g1]->ant = VP_direto[g2];
          VP_direto[g1]->prox  = VP_direto[PG2]; VP_direto[PG2]->ant = VP_direto[g1];

          //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
          VP_ultOpt = VP_geneIni->prox->id;

          recalcDist ();
          return;
      }
      else if (VP_direto[g2]->prox->id==g1)
      {
         int AG2 = VP_direto[g2]->ant->id;
         int PG1 = VP_direto[g1]->prox->id;

         VP_direto[AG2]->prox = VP_direto[g1]; VP_direto[g1]->ant = VP_direto[AG2];
         VP_direto[g1]->prox  = VP_direto[g2]; VP_direto[g2]->ant = VP_direto[g1];
         VP_direto[g2]->prox  = VP_direto[PG1]; VP_direto[PG1]->ant = VP_direto[g2];

         //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         recalcDist ();
         return;
      }

      for (tempGene=prox(VP_geneIni); tempGene->id!=VP_geneIni->id; tempGene = prox(tempGene))
      {
         //Encontrei um gene para troca
         if ((tempGene->id==g1)||(tempGene->id==g2))
         {
            if (!troca)
            {
               if(tempGene->id == g1)
               {
                  guarda = VP_direto[g2]->prox->id;
                  tempGene->ant->prox = VP_direto[g2];
                  VP_direto[g2]->prox = tempGene->ant;
                  tempGene->ant=tempGene->prox;
                  tempGene = VP_direto[g2];
               }
               else if(tempGene->id == g2)
               {
                  guarda = VP_direto[g1]->prox->id;
                  tempGene->ant->prox = VP_direto[g1];
                  VP_direto[g1]->prox = tempGene->ant;
                  tempGene->ant=tempGene->prox;
                  tempGene = VP_direto[g1];
               }
            }
            else
            {
               tempGene->prox = VP_direto[guarda];
               VP_direto[guarda]->ant = tempGene;
            }
              
            troca = !troca;
         }
            
         if(troca)
         {
            vizinho = tempGene->prox->id;
            tempGene->prox = tempGene->ant;
            tempGene->ant = VP_direto[vizinho];
         }

         VP_dist += Mapa->get_distancia(tempGene->ant->ori, tempGene->dest);
      }

      //Se o indivíduo foi modificado, o 2opt ou 3 opt deve começar do início
      VP_ultOpt = VP_geneIni->prox->id;

      VP_dist += Mapa->get_distancia(VP_geneIni->ant->ori, VP_geneIni->dest);

      //checa((char *)"troca");
   }

   int *criaGene (int ini, int tam)
   {
      int *novoGene;
      TGene *g = VP_direto[ini];

      novoGene = (int *)malloc(tam*sizeof(int));

      for (int i=0; i<tam; i++)
      {
         novoGene[i] = g->id;
         g = g->prox;
      }
      return novoGene;
   }

   public:
      TIndividuo ()
      {
         VP_qtdeGeneAlloc = 0;
      }

      //As alocações de memória não são realizadas no construtor
      //Por isso a necessidade do contador e da função
      void destroi()
      {
         for (int i=0; i<VP_qtdeGeneAlloc; i++) free (VP_direto[i]);
         if (VP_qtdeGeneAlloc>0)                free (VP_direto);  
      }

      TMapaGenes *Mapa;

      TGene *get_ini () { return VP_geneIni; }

      TGene *prox (TGene *gene) { return gene->prox; }
      TGene *prox (int idG) { return VP_direto[idG]->prox; }

      TGene *ant (TGene *gene) { return gene->ant; }
      TGene *ant (int idG) { return VP_direto[idG]->ant; }

      int get_qtdeGenes() { return VP_qtdeGenes; }
      
      bool get_is2Opt () { return VP_ultOpt==0; }
      void set_is2Opt (bool is2Opt) { VP_ultOpt = (is2Opt)?0:VP_geneIni->prox->id; }

      void checa (char *label)
      {
         TGene *tempG;
         printf ("\n%s\npor prox : ", label);
         for (tempG=VP_direto[0]; tempG->prox->id!=0; tempG = tempG->prox) printf ("%d - ", tempG->id);
         printf ("%d - %d\n", tempG->id, tempG->prox->id);

         printf ("\npor ant : ");
         for (tempG=VP_direto[0]; tempG->ant->id!=0; tempG = tempG->ant) printf ("%d - ", tempG->id);
         printf ("%d - %d\n\n", tempG->id, tempG->ant->id);
      }

      //Cria um novo indivíduo
      void novo ()
      {
         int i;

         VP_qtdeGenes = Mapa->get_qtdeGenes();
         
         VP_direto = (TGene **) malloc(VP_qtdeGenes*sizeof(TGene *));

         VP_dist = 0;

         VP_geneIni = (TGene *) malloc(sizeof(TGene));
         VP_qtdeGeneAlloc++;

         VP_geneIni->id   = 0;
         VP_geneIni->ori  = 0;
         VP_geneIni->dest = 0;
         
         VP_direto[0] = VP_geneIni;

         for (i=1; i<VP_qtdeGenes; i++)
         {
            VP_direto[i] = (TGene *) malloc(sizeof(TGene));
            VP_qtdeGeneAlloc++;
   
            VP_direto[i]->id = i;
            
            //No primero nível da recursividade um gene representa apenas uma cidade,
            //desta forma, tanto a origem como o destino é o próprio gene
            VP_direto[i]->ori = i; 
            VP_direto[i]->dest = i; 
            VP_direto[i-1]->prox = VP_direto[i];
            VP_direto[i]->ant = VP_direto[i-1];
           
            VP_dist += Mapa->get_distancia(VP_direto[i-1]->ori, VP_direto[i]->dest);
         }

         //Fechando o ciclo
         VP_direto[VP_qtdeGenes-1]->prox = VP_direto[0];
         VP_direto[0]->ant = VP_direto[VP_qtdeGenes-1];
           
         //Se o indivíduo foi modificado, o 2opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         VP_dist += Mapa->get_distancia(VP_direto[VP_qtdeGenes-1]->ori, VP_direto[0]->dest);

         //checa((char *)"novo");
      }


      //Cria um novo indivíduo
      void novo (int qtde, TTipoConvercao **tabConvercao)
      {
         int i;

         VP_qtdeGenes = qtde;
         
         VP_direto = (TGene **) malloc(VP_qtdeGenes*sizeof(TGene *));

         VP_dist = 0;

         VP_geneIni = (TGene *) malloc(sizeof(TGene));
         VP_qtdeGeneAlloc++;

         VP_geneIni->id   = tabConvercao[0]->id;
         VP_geneIni->ori  = tabConvercao[0]->ori;
         VP_geneIni->dest = tabConvercao[0]->dest;
         
         VP_direto[0] = VP_geneIni;

         for (i=1; i<VP_qtdeGenes; i++)
         {
            VP_direto[i] = (TGene *) malloc(sizeof(TGene));
            VP_qtdeGeneAlloc++;
   
            VP_direto[i]->id   = tabConvercao[i]->id;
            
            //No primero nível da recursividade um gene representa apenas uma cidade,
            //desta forma, tanto a origem como o destino é o próprio gene
            VP_direto[i]->ori  = tabConvercao[i]->ori;
            VP_direto[i]->dest = tabConvercao[i]->dest;

            VP_direto[i-1]->prox = VP_direto[i];
            VP_direto[i]->ant = VP_direto[i-1];
           
            VP_dist += Mapa->get_distancia(VP_direto[i-1]->ori, VP_direto[i]->dest);
         }

         //Fechando o ciclo
         VP_direto[VP_qtdeGenes-1]->prox = VP_direto[0];
         VP_direto[0]->ant = VP_direto[VP_qtdeGenes-1];
           
         //Se o indivíduo foi modificado, o 2opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         VP_dist += Mapa->get_distancia(VP_direto[VP_qtdeGenes-1]->ori, VP_direto[0]->dest);

         //checa((char *)"novo");
      }

      //distância total
      double get_distancia () { return VP_dist; }

      char *toString ()
      {

         TGene *tempGene;
         char *result;
         char temp[20];

         sprintf(temp, "%d", VP_geneIni->id);
         result = (char *)malloc(7*VP_qtdeGenes*sizeof(char));
         strcpy (result, temp);

         for (tempGene=prox(VP_geneIni); tempGene->id!=VP_geneIni->id; tempGene = prox(tempGene))
         {
            sprintf (temp, " - %d", tempGene->id);
            strcat (result, temp);            
         }
         return result;
      }

      char *toString (int init)
      {

         TGene *tempGene;
         char *result;
         char temp[20];

         sprintf(temp, "%d", VP_geneIni->id+init);
         result = (char *)malloc(7*VP_qtdeGenes*sizeof(char));
         strcpy (result, temp);

         for (tempGene=prox(VP_geneIni); tempGene->id!=VP_geneIni->id; tempGene = prox(tempGene))
         {
            sprintf (temp, " - %d", tempGene->id+init);
            strcat (result, temp);            
         }
         return result;
      }      

      //Embaralha os genes de um determinado individuo
      void embaralha ()
      { 
         int i;
         int rd1;
         int rd2;

         int G1A;
         int G1P;
         
         int G2A;
         int G2P;                  

         for (i=0; i<(VP_qtdeGenes/2); i++)
         {
            rd1 = (rand()%(VP_qtdeGenes-1))+1;
            G1A = VP_direto[rd1]->ant->id;
            G1P = VP_direto[rd1]->prox->id;

            rd2 = (rand()%(VP_qtdeGenes-1))+1;
            G2A = VP_direto[rd2]->ant->id;
            G2P = VP_direto[rd2]->prox->id; 

            //Não posso fazer o anterior de rd1 = anterior de rd2, pois ficaria anterior de rd1 = rd1
            if (VP_direto[rd1]->prox->id == VP_direto[rd2]->id)
            {
               VP_direto[rd1]->prox = VP_direto[G2P]; VP_direto[G2P]->ant = VP_direto[rd1];
               VP_direto[rd1]->ant = VP_direto[rd2]; VP_direto[rd2]->prox = VP_direto[rd1];
               VP_direto[rd2]->ant = VP_direto[G1A]; VP_direto[G1A]->prox = VP_direto[rd2];
            }
            //Não posso fazer o posterior de rd1 = posterior de rd2, pois ficaria posterior de rd1 = rd1
            else if (VP_direto[rd1]->ant->id == VP_direto[rd2]->id)
            {
               VP_direto[rd1]->ant = VP_direto[G2A]; VP_direto[G2A]->prox = VP_direto[rd1];
               VP_direto[rd1]->prox = VP_direto[rd2]; VP_direto[rd2]->ant = VP_direto[rd1];
               VP_direto[rd2]->prox = VP_direto[G1P]; VP_direto[G1P]->ant = VP_direto[rd2];
            }
            else
            {
               VP_direto[rd1]->prox = VP_direto[G2P]; VP_direto[G2P]->ant = VP_direto[rd1];
               VP_direto[rd1]->ant = VP_direto[G2A]; VP_direto[G2A]->prox = VP_direto[rd1];
               VP_direto[rd2]->prox = VP_direto[G1P]; VP_direto[G1P]->ant = VP_direto[rd2];
               VP_direto[rd2]->ant = VP_direto[G1A]; VP_direto[G1A]->prox = VP_direto[rd2];
            }
         }
        
         //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         recalcDist ();
         
         //checa((char *)"embaralha");
      } 

      //copia a sequência de genes do indivíduo para um array 
      void toArray (int *indArray)
      {
	 int i;
         TGene *tempG;

         for (i=0, tempG = VP_geneIni; i<VP_qtdeGenes; i++, tempG =  prox(tempG)) 
	    indArray[i] = tempG->id;
      }

      //copia a sequêcia de um array para o indivíduo
      void getArray (int *indArray)
      {
 	 int i;
	 VP_dist = 0;
	 VP_geneIni = VP_direto[indArray[0]];

	 for (i=1; i<VP_qtdeGenes; i++)
         {
            VP_direto[indArray[i-1]]->prox = VP_direto[indArray[i]];
            VP_direto[indArray[i]]->ant = VP_direto[indArray[i-1]];
           
            VP_dist += Mapa->get_distancia(VP_direto[indArray[i-1]]->ori, VP_direto[indArray[i]]->dest);
         }

         //Fechando o ciclo
         VP_direto[indArray[VP_qtdeGenes-1]]->prox = VP_direto[indArray[0]];
         VP_direto[indArray[0]]->ant = VP_direto[indArray[VP_qtdeGenes-1]];
           
         //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;

         VP_dist += Mapa->get_distancia(VP_direto[indArray[VP_qtdeGenes-1]]->ori, VP_direto[indArray[0]]->dest);

         //checa((char *)"getArray");

      }

      //Obtém o melhor indivíduo possível, pela realização de todas as combinações
      void melhorPossivel ()
      {
         int *indTemp;
         double melhorDist;
	 
         indTemp = (int *)malloc(VP_qtdeGenes*sizeof(int));
	 melhorDist = VP_dist;

	 toArray(indTemp);

         //A combinação é realizada de forma recursiva
         //pelo método combina
         combina (VP_geneIni->prox, indTemp, &melhorDist);

         if (melhorDist<VP_dist) getArray(indTemp);  
         free(indTemp);   

         //checa((char *)"melhorPossivel");
      }

      void debug ()
      {
         int i;
         for (i=0; i<VP_qtdeGenes; i++)
            printf ("Gene: %d, prox: %d, ant: %d\n", VP_direto[i]->id, VP_direto[i]->prox->id, VP_direto[i]->ant->id);
      }

      void recalcDist ()
      {
         TGene *tempGene;
         VP_dist = 0;

         for (tempGene=prox(VP_geneIni); tempGene->id!=VP_geneIni->id; tempGene = prox(tempGene))
             VP_dist += Mapa->get_distancia(ant(tempGene)->ori, tempGene->dest);
 
         VP_dist += Mapa->get_distancia(ant(VP_geneIni)->ori, VP_geneIni->dest);
      }

      //Realiza o 2-Opt
      void mutacao (int profundidade)
      {
         TGene *g1;
         TGene *g2;
         double tot;
      //   printf ("vou fazer 2opt ? %d\n", profundidade);
         if (VP_ultOpt == 0) return;
         VP_ultOpt = 0;
//printf("\nfiz 2opt profund. %d\n", profundidade);
         for (g1= VP_geneIni->prox; (g1->prox->prox->id!=0)&&(g1->id!=0); g1=g1->prox)
            for (g2=g1->prox; g2->prox->id!=0; g2=g2->prox)
            {  
//printf("for g2 %d - %d\n", g1->id, g2->id);
               tot = calcDistTotTroca (g1->id, g2->id);

               if(tot<VP_dist) 
               {
                  troca(g1->id, g2->id);
                  g1 =  VP_direto[g2->id];
//printf("voltando ao inico");
                  VP_ultOpt = g1->id;
               }
            }
            
           // if (VP_ultOpt>0) mutacao (profundidade);
//printf ("sai do 2opt");

         //Se saiu do for, é porque não tem o que o 2opt melhorar mais
         //VP_ultOpt = 0;

         //checa((char *)"mutacao");
      }

      TGene *get_gene (int id) { return VP_direto[id]; }

      //Cross Over, transforma o indivíduo em um descendente de outros 2
      void descendente (TIndividuo parceiro1, TIndividuo parceiro2)
      {
         bool *controle;
         TGene *gPar1;
         TGene *gPar2;

         TGene *g1;
         TGene *g2;

         int pivo = (rand()%(VP_qtdeGenes-1))+1;

         //Servirá de controle para os genes já utilzados
         controle = (bool *)calloc(VP_qtdeGenes, sizeof(bool));
         
         gPar1 = parceiro1.get_gene(pivo);
         gPar2 = parceiro2.get_gene(pivo);

         controle[pivo]=true;

         g1 = VP_direto[pivo];
         g2 = VP_direto[pivo];

         bool esq=true;
         bool dir=true;

         int i=1;

         //Enquanto não posicionar todos os genes
         while (i<VP_qtdeGenes)
         {
            if (esq)
            {
               if (gPar1->id==0) esq = false;
               else
               {
                  gPar1 = parceiro1.ant(gPar1->id);
                  if (!controle[gPar1->id])
                  {
                     g1->ant = VP_direto[gPar1->id];
                     VP_direto[gPar1->id]->prox = g1;
                     g1 = VP_direto[gPar1->id];
                     controle[gPar1->id] = true;
                     i++;
                  }
               }
            }

            if (dir)
            {
               if (gPar2->prox->id==0) dir = false;
               else
               {
                  gPar2 = parceiro2.prox(gPar2->id);
                  if (!controle[gPar2->id])
                  {
                     g2->prox = VP_direto[gPar2->id];
                     VP_direto[gPar2->id]->ant = g2;
                     g2 = VP_direto[gPar2->id];
                     controle[gPar2->id] = true;
                     i++;
                  }
               }
            }
            
            if ((!esq)&&(esq==dir))
            {
               gPar1 = parceiro1.ant(gPar1->id);
               if (!controle[gPar1->id])
               {
                  g1->ant = VP_direto[gPar1->id];
                  VP_direto[gPar1->id]->prox = g1;
                  g1 = VP_direto[gPar1->id];
                  controle[gPar1->id] = true;
                  i++;
               }

               gPar2 = parceiro2.prox(gPar2->id);
               if (!controle[gPar2->id])
               {
                  g2->prox = VP_direto[gPar2->id];
                  VP_direto[gPar2->id]->ant = g2;
                  g2 = VP_direto[gPar2->id];
                  controle[gPar2->id] = true;
                  i++;
               }
            }
         }
         g1->ant = g2;
         g2->prox = g1;
         free(controle);

         //Se o indivíduo foi modificado, o 2opt ou 3opt deve começar do início
         VP_ultOpt = VP_geneIni->prox->id;
         recalcDist ();

         //checa((char *)"descendente");
      }

      TTipoConvercao **criaTabelaConversao(TIndividuo *melhor, int *qtde)
      {
         TTipoConvercao **tabConvercao;
         int *novoGene;
         TTipoConvercao *tabTemp;
	 TTipoConvercao *tab;
	 TTipoConvercao *tabIni;

         int tamGene=1;
         TGene *g;
         TGene *tmp;

         int gIni;                   

         *qtde=0;

         gIni = 0;

         for (g=melhor->get_gene(0); (g->id!=0)||(*qtde==0); g=g->prox)
         {
            tmp = VP_direto[g->id];

            if ((g->prox->id == tmp->prox->id)&&(g->prox->id!=0)) tamGene++;
            else
            {
               novoGene = criaGene (gIni, tamGene); 

               tabTemp = (TTipoConvercao *) malloc(sizeof(TTipoConvercao));
               tabTemp->gene = novoGene; 
	       tabTemp->id = *qtde;
               tabTemp->ori = g->ori;
               tabTemp->dest = VP_direto[gIni]->dest;
	       tabTemp->qtde = tamGene;

               if(*qtde==0) tabIni = tabTemp;
               else
               {
                  tab->prox = tabTemp;
                  tabTemp->ant = tab;
               }
               tab = tabTemp;
               gIni = g->prox->id;
	       tamGene = 1;
               *qtde=*qtde+1;
            }
         }

         tabConvercao = (TTipoConvercao **)malloc(*qtde*sizeof(TTipoConvercao *));

         int i;

         for(i=0, tabTemp=tabIni; i<*qtde; i++, tabTemp=tabTemp->prox)
            tabConvercao[i] = tabTemp;

         return tabConvercao;
      }

      void converte(TIndividuo *geneConv, TTipoConvercao **tabConvercao)
      {
         int id;
         TGene *gA;
         TGene *gTemp;
         int qtde = 0;

         gA = VP_direto[0];

         for (gTemp=geneConv->get_gene(0); (gTemp->id!=0)||(qtde==0); gTemp=geneConv->prox(gTemp->id))
         {
            qtde = tabConvercao[gTemp->id]->qtde;
            for (int i=0; i<qtde; i++)
            {
               id = tabConvercao[gTemp->id]->gene[i];
               if (id==0) continue;
               gA->prox = VP_direto[id];
               VP_direto[id]->ant = gA;
	       gA = gA->prox;
            }
	 }

         gA->prox = VP_direto[0];
         VP_direto[0]->ant = gA;
         recalcDist ();

         //checa((char *)"converte");
      }

};

/*********************************************************
classe de População, conterá todos os indivíduos e gerará
as novas gerações
**********************************************************/
class TPopulacao
{
   private:
      double VP_somaDistancias;
      int VP_tamanho;
      
      int VP_maxGeracao;
      int VP_printParcial;
      int VP_percentMutacao;
      int VP_percentManipulacao;
      int VP_percentReducao;

      int VP_profundidade;

      TIndividuo **VP_individuos;
      TMapaGenes *VP_mapa;

      void liberaTabConversao(TTipoConvercao **tabConvercao, int qtdeGenes)
      {
         for(int i=qtdeGenes-1; i>=0; i--)
         {
            free (tabConvercao[i]->gene);
            free (tabConvercao[i]);
         }

         free (tabConvercao);
      }

   public:
      //***********************
      //   Propriedades
      //***********************
      void set_maxGeracao (int val) { VP_maxGeracao = val; }
      int get_maxGeracao () { return VP_maxGeracao; }

      void set_printParcial (int val) { VP_printParcial = val; }
      int get_printParcial () { return VP_printParcial; }

      void set_percentMutacao (int val) { VP_percentMutacao = val; }
      int get_percentMutacao () { return VP_percentMutacao; }

      void set_percentManipulacao (int val) { VP_percentManipulacao = val; }
      int get_percentManipulacao () { return VP_percentManipulacao; }

      void set_percentReducao (int val) { VP_percentReducao = val; }
      int get_percentReducao () { return VP_percentReducao; }

      //*******************************
      //   Construtor e Destrutor
      //*******************************
      TPopulacao (int tamanho, TMapaGenes *mapa, int profundidade = 0)
      {
         if (profundidade>0)
         {
            if (tamanho<20) tamanho = 20;
            if (tamanho>2*mapa->get_qtdeGenes()) tamanho = 2*mapa->get_qtdeGenes();
         }

	 VP_profundidade = profundidade;
         VP_somaDistancias = 0;

         VP_mapa = mapa;

         VP_tamanho = tamanho;

         VP_individuos = (TIndividuo **)malloc(VP_tamanho*sizeof(TIndividuo *));         
         for (int i=0; i<VP_tamanho; i++)
         {
            VP_individuos[i] = new TIndividuo();
            VP_individuos[i]->Mapa = VP_mapa;
            VP_individuos[i]->novo();

            //Adicionado, para que o melhor indivíduo seja o mesmo
            //Melhor indivíduo anterior
            if (i>0) VP_individuos[i]->embaralha();

            VP_somaDistancias += VP_individuos[i]->get_distancia();
         }         
      }

      TPopulacao (bool is2Opt, int tamanho, TMapaGenes *mapa, int qtde, TTipoConvercao **tabConvercao, int profundidade = 0)
      {
         //printf ("is2opt = %d\n", is2Opt);
         if (profundidade>0)
         {
            if (tamanho<20) tamanho = 20;
            if (tamanho>2*qtde) tamanho = 2*qtde;
         }

	 VP_profundidade = profundidade;
         VP_somaDistancias = 0;

         VP_mapa = mapa;

         VP_tamanho = tamanho;

         VP_individuos = (TIndividuo **)malloc(VP_tamanho*sizeof(TIndividuo *));         
         for (int i=0; i<VP_tamanho; i++)
         {
            VP_individuos[i] = new TIndividuo();
            VP_individuos[i]->Mapa = VP_mapa;
            VP_individuos[i]->novo(qtde, tabConvercao);

            //Adicionado, para que o melhor indivíduo seja o mesmo
            //Melhor indivíduo anterior
            if (i>0) VP_individuos[i]->embaralha();
            else VP_individuos[i]->set_is2Opt(is2Opt);

            VP_somaDistancias += VP_individuos[i]->get_distancia();
         }         
      }

      ~TPopulacao()
      {
         for (int i=VP_tamanho-1; i>=0; i--)
         {
            VP_individuos[i]->destroi();
            free (VP_individuos[i]);
         }
         free (VP_individuos);
      }

      //***********************
      //   Métodos
      //***********************
      char *toString ()
      {
         char *result;
         char temp[25];
         char *str;

         result = (char *)malloc(VP_tamanho*(8*VP_individuos[0]->get_qtdeGenes()+25)*sizeof(char));
         str = VP_individuos[0]->toString();
         strcpy (result, str);
         free (str);
         sprintf (temp, " (%f)", VP_individuos[0]->get_distancia());
         strcat (result, temp);           

         for (int i=1; i<VP_tamanho; i++)
         {
            strcat (result, "\n"); 
            str = VP_individuos[i]->toString();
            strcat (result, str);
            free (str);   
            sprintf (temp, " (%f)", VP_individuos[i]->get_distancia());
            strcat (result, temp);           
         }
         return result;
      }     

      void ordena () 
      { 
         qsort (VP_individuos, VP_tamanho, sizeof(TIndividuo *), compare);
      }  

      void novaGeracao()
      {
         int qtdePercent = VP_tamanho * VP_percentManipulacao/100;
         int qtdeAlt = 0;
         double difDist;

         int anterior = 0;
         //Substituindo todos que tem diferença de distância 0 para o anterior
         for (int i=1; i<VP_tamanho; i++)
         {
            difDist = VP_individuos[i]->get_distancia() - VP_individuos[anterior]->get_distancia();
            if(difDist == 0)
            {
               int parceiro1=i;
               int parceiro2=i;

               while (parceiro1==i) parceiro1 = (rand()%VP_tamanho); //Não pode fazer parte do crossover

               // Não permite que o indivíduo cruze com individuo igual a ele
               while ((parceiro2==i)||(parceiro2==parceiro1)) parceiro2 = (rand()%VP_tamanho);

               VP_somaDistancias -= VP_individuos[i]->get_distancia();

               if (VP_individuos[parceiro1]->get_distancia()==VP_individuos[parceiro2]->get_distancia())
                  VP_individuos[i]->embaralha();
               else
                  VP_individuos[i]-> descendente(*VP_individuos[parceiro1],  *VP_individuos[parceiro2]);

               VP_somaDistancias += VP_individuos[i]->get_distancia();

	       qtdeAlt++;
               if (qtdeAlt>=qtdePercent) break;
            }
            else
               anterior = i;           
         }

         //Completa o crossover com os piores.
         //Obs.: um dos que tem distância 0 pode tabém ser um dos piores, mas, 
         //      não é preciso se preocupar com isso
         for (int i=VP_tamanho-qtdePercent+qtdeAlt; i<VP_tamanho; i++)
         {
            int parceiro1=i;
            int parceiro2=i;

            while (parceiro1==i) parceiro1 = (rand()%VP_tamanho); //Não pode fazer parte do crossover

            // Não permite que o indivíduo cruze com individuo igual a ele
            while ((parceiro2==i)||(parceiro2==parceiro1)) parceiro2 = (rand()%VP_tamanho);

            VP_somaDistancias -= VP_individuos[i]->get_distancia();
               
            if (VP_individuos[parceiro1]->get_distancia()==VP_individuos[parceiro2]->get_distancia())
               VP_individuos[i]->embaralha();
            else
               VP_individuos[i]-> descendente(*VP_individuos[parceiro1],  *VP_individuos[parceiro2]);

            VP_somaDistancias += VP_individuos[i]->get_distancia();
         }
      } 

      void mutacaoAGRecursivo(int escolha)
      {
         //Tranosformará as sequências, do indivíduo, que coincidem com o melhor em um gene
         TTipoConvercao **tabConvercao;
         int qtdeGenes;
         double esforco;
         double fatorial;
	 TMapaGenes tempMapa;

         TIndividuo *melhor = VP_individuos[0];
         TIndividuo *manipulado = VP_individuos[escolha];

	 //Verifico se o esforço
         //(poulação * % percent de manipulação)* Número de gerações * (população * %mutação)
	 esforco = (VP_tamanho * VP_percentMutacao/100) * VP_maxGeracao * (VP_tamanho * VP_percentManipulacao/100);

         //So o fatorial for <= 0 é porque houve overflow
         fatorial = fatorialStirling(manipulado->get_qtdeGenes());

         //é maior do que todas as combinações (fatorial) da quantidade de gene
	 //Se a combinação form maior, realiza o AG com os novos genes.
         if ((fatorial>0)&&(esforco>=fatorial))
         {
            
            //Quando o esforço é maior que as combinações
            //O melhor indivuo por combinação
	    manipulado->melhorPossivel();
            
            if (manipulado->get_distancia() < melhor->get_distancia())
               printf ("Melhorou na combinação\n");

            return;
         }

	 // Novos Genes para AG
         //********************
	 // tentará reduzir a quantidade de gene.
         tabConvercao = manipulado->criaTabelaConversao(melhor, &qtdeGenes);
         // Se não reduzir, não faz nada
         //Obs.: qtdeGenes<=2, pois se for 1, indivíduo é igual ao melhor
         //      se for 2, não posso trocar as posições, já que sempre inicia
         //      pelo gene 0. Logo, posso sair sem fazer nada

         if ((qtdeGenes >= VP_percentReducao*0.01*manipulado->get_qtdeGenes())||(qtdeGenes<=2))
         {
            //Libera memória
            liberaTabConversao(tabConvercao, qtdeGenes);
	    return;
         }

         //Reduzindo o tamanho da população a cada integração
         //float per = (float) qtdeGenes/manipulado->get_qtdeGenes();
         //int newTam = VP_tamanho*per;

         //realizar o AG
         // melhor->get_is2Opt() : se o melhor indivíduo já foi melhorado a exaustão pelo 2opt
         // Não hpa porquer tentar melhor sua redução

         TPopulacao populacao (melhor->get_is2Opt(), VP_tamanho, VP_mapa, qtdeGenes, tabConvercao, VP_profundidade+1);

         //Configurando a população para execução do AG
	 populacao.set_maxGeracao (VP_maxGeracao);
	 populacao.set_printParcial (VP_printParcial);
         populacao.set_percentMutacao (VP_percentMutacao);
         populacao.set_percentManipulacao (VP_percentManipulacao);
         populacao.set_percentReducao (VP_percentReducao);

         //A execução do algoritmo genético
         populacao.algaritmoGenetico ();

         //converte o melhor indivíduo do AG para o indivíduo atual
         manipulado->converte(populacao.get_melhor(), tabConvercao);

         //Libera memória
         liberaTabConversao(tabConvercao, qtdeGenes);

      }

      void mutacao()
      {
         int qtdePercent = VP_tamanho * VP_percentMutacao/100;

	//A mutação será:
        //O primero indivíduo sofre mutação 2-opt
        VP_somaDistancias -=  VP_individuos[0]->get_distancia();
        VP_individuos[0]->mutacao(VP_profundidade);
        VP_somaDistancias +=  VP_individuos[0]->get_distancia();

        //Os demais, aleatóriamente sofrerão uma mutação de melhora recursiva
        for (int i=1; i<qtdePercent; i++)
        {
            int escolha = (rand()%(VP_tamanho-1)+1);
            VP_somaDistancias -= VP_individuos[escolha]->get_distancia();
            mutacaoAGRecursivo(escolha);
            VP_somaDistancias += VP_individuos[escolha]->get_distancia();

            if (VP_individuos[escolha]->get_distancia()<VP_individuos[0]->get_distancia())
            {
               printf("\nmelhorou - profundidade %d\n", VP_profundidade+1);
               TIndividuo melhor = *VP_individuos[0];
	       *VP_individuos[0] = *VP_individuos[escolha];
               *VP_individuos[escolha] = melhor;
               //Se melhorou, é necessário voltar a fazer o 2opt
               VP_individuos[0]->set_is2Opt (false); 
            }
        }
      } 

      void embaralha()
      {
         for (int i = 1; i<VP_tamanho; i++)
            VP_individuos[i]->embaralha();
      }

      TIndividuo *get_melhor() { return VP_individuos[0]; }
      TIndividuo *get_pior() { return VP_individuos[VP_tamanho-1]; }
      TIndividuo *get_individuo(int index) { return VP_individuos[index]; }
      double distanciaMedia () { return VP_somaDistancias/VP_tamanho; }

      //******************************************************
      //   O    A L G O R I T M O   G E N É T I C O 
      //******************************************************
      void algaritmoGenetico ()
      {
         TIndividuo *melhor;
         TIndividuo *pior;

         double m = infinito;
         ordena();

         for(int i=0; i<VP_maxGeracao; i++)
         {
            novaGeracao ();
            mutacao ();
            ordena();
            
	    melhor = get_melhor();
            if((VP_printParcial)&&(VP_profundidade==0))
            {
               time(&sysTime2);
               pior = get_pior();
               printf ("Geração %d. Melhor : %f - Pior: %f - Media: %f (%f)\n",i, melhor->get_distancia(), pior->get_distancia(), distanciaMedia(), difftime(sysTime2, sysTime1));
            }

            if (m<=melhor->get_distancia())
            {
               if (VP_profundidade>0) break;               
            } 
            
            m = melhor->get_distancia();
         }

         if((VP_printParcial)&&(VP_profundidade==0))
         {
            //Imprimindo o melhor indivíduo    
            char *str = melhor->toString();
            printf ("%s (%f)\n", str, melhor->get_distancia());
            free (str);
         }

      }
};

int compare(const void *x, const void *y)
{
  return ((*(TIndividuo **)x)->get_distancia() - (*(TIndividuo **)y)->get_distancia());
}

/*********************************************************
Classe de confguração, utilizada para acelerar o processo
de testes, estava ruim compilar para cada configuração e
estava ruim digitar todos os parâmetros de confguração
**********************************************************/
class TConfig
{
   private:

   void leInfo(xmlDocPtr doc, xmlNode * a_node)
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
         }

         leInfo(doc, cur_node->children);
      }
   }

   public:
      /*******************************************
           Parâmetros de teste do sistema
      ********************************************/
      int tamPopulacao;
      int maxGeracao;
      int percentManipulacao;   //percentual de manipulação do indivíduo. (exclusão / cruzamento)
      int percentMutacao;       //percentual de mutação
      int printParcial;         //Se inprime informações intermediárias
      int percentReducao;       //percentual de redução do gene a cada recursvidade

      /*******************************************************
           Os valores padrões são os utilizados no artigo
      ********************************************************/
      TConfig()
      {
         tamPopulacao = 200;
         maxGeracao = 300;
         percentManipulacao = 30;
         percentMutacao = 20;
         percentReducao = 50;
      }

      void carregaDoArquivo(char *nomeArquivo)
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
};

/*************************************************************************
 Principal, o arquivo xml a ser carregado deve ser passado como parâmetro

 Linha de execução
 1) executa o programa com as confgurações padrões
    tsp 

 2) executa o programa com interação sobre os valores de configuração
    tsp -c  

 3) executa o programa lendo um arquivo de configuração em xml
    tsp -c <arquivo de confguração> 
**************************************************************************/
int main(int argc, char **argv)
{
   TMapaGenes mapa ;
   TIndividuo *melhor;
   TConfig config; 

   if (argc < 2)
      return(1);

   LIBXML_TEST_VERSION

   mapa.carregaDoArquivo (argv[1]);

   if (argc>2)
   {
      if ((argv[2][0] == '-')&&(argv[2][1] == 'c'))
      {
         if (argc > 3)
         {
            config.carregaDoArquivo(argv[3]);
         }
         else
         {
            printf("Mostra informações parciais: "); scanf("%d", &config.printParcial);
            printf("Tamanho da população: "); scanf("%d", &config.tamPopulacao);
            printf("Número máximo de gerações: "); scanf("%d", &config.maxGeracao);
            printf("Percentual de manipuação da população: "); scanf("%d", &config.percentManipulacao);
            printf("Percentual de mutação: "); scanf("%d", &config.percentMutacao);
            printf("Percentual de Reducao: "); scanf("%d", &config.percentReducao);
         }
      }
      else
         return 1;
   }

    time(&sysTime1);

    //Criando a população
    TPopulacao populacao (config.tamPopulacao, &mapa);

    //Configurando a população para execução do AG
    populacao.set_maxGeracao (config.maxGeracao);
    populacao.set_printParcial (config.printParcial);
    populacao.set_percentMutacao (config.percentMutacao);
    populacao.set_percentManipulacao (config.percentManipulacao);
    populacao.set_percentReducao (config.percentReducao);

    //A execução do algoritmo genético
    populacao.algaritmoGenetico ();

    time(&sysTime2);
    printf ("Temmpo de execução %f\n", difftime(sysTime2, sysTime1));

    //Imprimindo o melhor indivíduo    
    melhor = populacao.get_melhor();
    char *str = melhor->toString(1);
    printf ("%s (%f)\n", str, melhor->get_distancia());
    free (str);

    str = populacao.get_individuo(1)->toString(1);
    printf ("\n\n1: %s (%f)\n", str, populacao.get_individuo(1)->get_distancia());
    free (str);

    str = populacao.get_individuo(2)->toString(1);
    printf ("\n\n2: %s (%f)\n", str, populacao.get_individuo(2)->get_distancia());
    free (str);

    str = populacao.get_individuo(3)->toString(1);
    printf ("\n\n3: %s (%f)\n", str, populacao.get_individuo(3)->get_distancia());
    free(str);

    str = populacao.get_pior()->toString(1);
    printf ("\n\n pior: %s (%f)\n", str, populacao.get_pior()->get_distancia());
    free(str);
    return 0;
}
#else
int main(void) {
    fprintf(stderr, "Suporte a árvore xml não compilado\n");
    exit(1);
}
#endif
