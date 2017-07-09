/*
*  main.cpp
*
*  Módulo pricipal do estudo TSP
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 05/05/2017
*  Datas de Modificações:
*
* Programa de Algorítmo Genetico para resolver tsp de forma recursiva
* Necessário instalação da biblioteca libxml
* Instalação da biblioteca no debian
*           $> sudo apt-get install libxml
*
*  Os parametros de entrada estão explicados no arquivo tsphelp.txt
*  A compilação pode ser realizada pelo comando make
*
*  Se não for definido parametros de entrada, serão utilizados os
*  parametros padrões, conforme abaixo.
*
*  número de execuções:
*  mutação:
*  cruzamento:
*  população:
*  % de mutação:
*  % de cruzamento:
*  recursividade:
*  % de melhora recursiva:
*  arquivo de saida:
*
*/
#include <iostream>
#include <string>
#include "config.hpp"
#include "ag.hpp"
#include "tsp.hpp"

using namespace std;

#ifdef LIBXML_TREE_ENABLED

int main(int argc, char *argv[])
{
   string nomeArqSaida;
   string cabecalho;
   TAlgGenetico *ag;
   TMapaGenes *mapa = new TMapaGenes();
   TConfig *config  = new TConfig();
   TArqLog *arqSaida;

   time_t tempo;
   struct tm *tlocal;
   char data[128];

   //parâmetros obrigatóros como entrada
   if (argc < 4)
   {
      cout << "Parâmetros obrigatóros:" << endl;
      cout << "\t 1 - Arquivo de instância TSP, no formato XML" << endl;
      cout << "\t 2 - Arquivo de configuração, no formato XML" << endl;
      cout << "\t 3 - Nome do arquivo de saída, resultados" << endl;
      return 1;
   }

   LIBXML_TEST_VERSION

   cout << "Caregando arquivo de configuração: " << argv[2] << endl;
   config->carregaDoArquivo(argv[2]);
   cout << "Arquivo " << argv[2] << " carregado." << endl;

   cout << "Caregando instância " << argv[1] << endl;
   mapa->carregaDoArquivo (argv[1]);
   cout << "Instância " << argv[1] << " carregada." << endl;

   TUtils::initRnd ();

   for (int countExec=0; countExec<config->numExec; countExec++)
   {
      /* para complemento do nome do arquivo de saida */
      tempo = time(0);
      tlocal = localtime(&tempo);
      strftime(data, 128, "%d_%m_%y_%H_%M_%S", tlocal);
      nomeArqSaida = argv[3];
      nomeArqSaida += "_";
      nomeArqSaida += to_string(countExec);
      nomeArqSaida += "_";
      nomeArqSaida += data;
      nomeArqSaida += ".txt";

      if (config->printParcial)
      {
         cout << "Execução " << countExec+1 << " / " << config->numExec << " iniciada. "<< endl;
         cout << "Arquivo de saida de resultados: " << nomeArqSaida << endl;
      }

      cabecalho   = "Instância;";
      cabecalho  +=  argv[1];
      cabecalho  +=  "\n";

      cabecalho  += "Execução;";
      cabecalho  += to_string(countExec+1);
      cabecalho  += " / ";
      cabecalho  += to_string(config->numExec);
      cabecalho  += "\n";

      cabecalho  += "Tamanho da Populaçap;";
      cabecalho  += to_string(config->tamPopulacao);
      cabecalho  += "\n";

      cabecalho  += "Máximo de gerações;";
      cabecalho  += to_string(config->maxGeracao);
      cabecalho  += "\n";

      cabecalho  += "Mutação;";
      cabecalho  += to_string(config->mutacao);
      cabecalho  += "\n";

      cabecalho  += "Cruzamento;";
      cabecalho  += to_string(config->cruzamento);
      cabecalho  += "\n";

      cabecalho  += "% Manpulação;";
      cabecalho  += to_string(config->percentManipulacao);
      cabecalho  += "\n";

      cabecalho  += "% Mutação;";
      cabecalho  += to_string(config->percentMutacao);
      cabecalho  += "\n";

      cabecalho  += "Print Parcial;";
      cabecalho  += to_string(config->printParcial);
      cabecalho  += "\n";

      cabecalho  += "Ativa Recursivo;";
      cabecalho  += to_string(config->percentMutacaoRecursiva);
      cabecalho  += "\n";

      cabecalho  += "Percentual de Redução;";
      cabecalho  += to_string(config->percentReducao);
      cabecalho  += "\n";

      cabecalho  += "Profundidade Máxima;";
      cabecalho  += to_string(config->profundidadeMaxima);
      cabecalho  += "\n";

      cabecalho  += "Percentual de Eltismo;";
      cabecalho  += to_string(config->percentElitismo);
      cabecalho  += "\n";

      cabecalho  += "Seleção para Cruzamento;";
      cabecalho  += to_string(config->selecao);
      cabecalho  += "\n";

      cabecalho  += "Forma de Seleção par Mutação;";
      cabecalho  += to_string(config->selIndMutacao);
      cabecalho  += "\n";

      arqSaida = new TArqLog(cabecalho, nomeArqSaida);
      ag = new TAlgGenetico(mapa, arqSaida);
      ag->setMutacao(config->mutacao);
      ag->setCruzamento(config->cruzamento);
      ag->setTamPopulacao(config->tamPopulacao);
      ag->setPrintParcial(config->printParcial);
      ag->setMaxGeracao(config->maxGeracao);
      ag->setPercentElitismo(config->percentElitismo);
      ag->setPercentMutacao(config->percentMutacao);
	   ag->setProfundidadeMaxima(config->profundidadeMaxima);
      ag->setSelecao(config->selecao);
      ag->setSelIndMutacao(config->selIndMutacao);
      ag->setPercentMutacaoRecursiva(config->percentMutacaoRecursiva);
      ag->setPercentReducao(config->percentReducao);
      ag->exec();

      arqSaida->addLinha("");
      arqSaida->addLinha("");

      delete arqSaida;
      delete ag;
   }

   delete config;
   delete mapa;
   return 0;
}
#else
int main(void)
{
    cout << "Suporte a árvore xml não compilado" << endl;
    exit(1);
}
#endif
