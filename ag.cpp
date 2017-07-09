/*
*  ag.cpp
*
*  Módulo responsável pela execução do algoritmo genético
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 06/05/2017
*  Datas de Modificações:
*
*/

#include "ag.hpp"
#include <iostream>

//Apenas letura
TMapaGenes *TAlgGenetico::getMapa() { return VP_Mapa; }
TArqLog *TAlgGenetico::getArqLog()  { return VP_ArqSaida; }
int TAlgGenetico::get_profundidade(){ return VP_profundidade; }

//propriedades
/**************/

void TAlgGenetico::setTamPopulacao (int val)            { VP_tamPopulacao            = val; }
void TAlgGenetico::setMaxGeracao (int val)              { VP_maxGeracao              = val; }
void TAlgGenetico::setPrintParcial (int val)            { VP_printParcial            = val; }
void TAlgGenetico::setPercentMutacao (int val)          { VP_percentMutacao          = val; }
void TAlgGenetico::setPercentManipulacao (int val)      { VP_percentManipulacao      = val; }
void TAlgGenetico::setProfundidadeMaxima(int val)       { VP_profundidadeMaxima      = val; }
void TAlgGenetico::setMutacao (int val)                 { VP_mutacao                 = val; }
void TAlgGenetico::setCruzamento (int val)              { VP_cruzamento              = val; }
void TAlgGenetico::setPercentElitismo (int val)         { VP_percentElitismo         = val; }
void TAlgGenetico::setSelecao (int val)                 { VP_selecao                 = val; }
void TAlgGenetico::setSelIndMutacao (int val)           { VP_selIndMutacao           = val; }
void TAlgGenetico::setPercentReducao (int val)          { VP_percentReducao          = val; }
void TAlgGenetico::setPercentMutacaoRecursiva (int val) { VP_percentMutacaoRecursiva = val; }

int TAlgGenetico::getPercentReducao ()          { return VP_percentReducao; }
int TAlgGenetico::getPercentMutacaoRecursiva () { return VP_percentMutacaoRecursiva; }
int TAlgGenetico::getTamPopulacao ()            { return VP_tamPopulacao; }
int TAlgGenetico::getMaxGeracao ()              { return VP_maxGeracao; }
int TAlgGenetico::getPrintParcial ()            { return VP_printParcial; }
int TAlgGenetico::getPercentMutacao ()          { return VP_percentMutacao; }
int TAlgGenetico::getPercentManipulacao ()      { return VP_percentManipulacao; }
int TAlgGenetico::getProfundidadeMaxima()       { return VP_profundidadeMaxima; }
int TAlgGenetico::getMutacao ()                 { return VP_mutacao; }
int TAlgGenetico::getCruzamento ()              { return VP_cruzamento; }
int TAlgGenetico::getPercentElitismo ()         { return VP_percentElitismo; }
int TAlgGenetico::getSelecao ()                 { return VP_selecao; }
int TAlgGenetico::getSelIndMutacao ()           { return VP_selIndMutacao; }

unsigned TAlgGenetico::getEntradaRec() { return VP_Entr_Rec; }
unsigned TAlgGenetico::getExecRec()    { return VP_Exec_Rec; }
unsigned TAlgGenetico::getCombinaRec() { return VP_Comb_Rec; }


//Métodos
/**********/
void TAlgGenetico::setTime (time_t sTime) {	sysTime1 = sTime; }

TAlgGenetico::TAlgGenetico (TMapaGenes *mapa, TArqLog *arqSaida)
{
   VP_Mapa = mapa;
   VP_ArqSaida = arqSaida;
   VP_profundidade = 0;
	VP_Entr_Rec 	 = 0;
	VP_Exec_Rec     = 0;
	VP_Comb_Rec     = 0;
}

TAlgGenetico::TAlgGenetico (TMapaGenes *mapa, TArqLog *arqSaida, int profundidade)
{
   VP_Mapa = mapa;
   VP_ArqSaida = arqSaida;
   VP_profundidade = profundidade;
	VP_Entr_Rec     = 0;
	VP_Exec_Rec     = 0;
	VP_Comb_Rec     = 0;
}

TAlgGenetico::~TAlgGenetico ()
{

}

double TAlgGenetico::get_melhor_dst()
{
   return VP_melhor_dist;
}

int TAlgGenetico::get_melhor_gera()
{
   return VP_melhor_gera;
}

void TAlgGenetico::exec(TIndividuo *indBase, TTabConversao *tabConversao, int extraPadrao)
{  
   //Limitando a profundidade da recursivdade
   if (get_profundidade() >= getProfundidadeMaxima()) return;
   
   TIndividuo *melhor;
   mutacao = new TMutacao(getMapa(), getArqLog(), getMutacao());
   cruzamento = new TCruzamento(getMapa(), getArqLog(), getCruzamento());
   
   int tam = getTamPopulacao();
   if (tam>2*tabConversao->get_quantidade()) tam=2*tabConversao->get_quantidade();
   
   TPopulacao *populacao = new TPopulacao (tam, getMapa(), getArqLog());
   vector<TTipoConversao> tc = tabConversao->get_vetConv();
   populacao->povoa(tc, extraPadrao);
   melhor = populacao->get_melhor();
   populacao->ordena();
   	   
   double mdAtu = melhor->get_distancia();
   
   for (int i=1; i<=getMaxGeracao(); i++)
   {
      populacao = iteracao (populacao, i);
		melhor = populacao->get_melhor();
		if(melhor->get_distancia()<mdAtu) mdAtu = melhor->get_distancia();
		else if (get_profundidade() > 0) break;
   }

   delete mutacao;
   delete cruzamento;

   tabConversao->converte(indBase, populacao->get_melhor());
   
   tc.clear();
   delete populacao;
}
   
void TAlgGenetico::exec()
{
   TIndividuo *melhor;
   mutacao = new TMutacao(getMapa(), getArqLog(), getMutacao());
   cruzamento = new TCruzamento(getMapa(), getArqLog(), getCruzamento());

   VP_melhor_dist = infinito;
   TPopulacao *populacao = new TPopulacao (getTamPopulacao(), getMapa(), getArqLog());
   populacao->povoa();
/*
//Descomente para teste de cruzamento
TIndividuo *i1 = populacao->get_individuo(0);
TIndividuo *i2 = populacao->get_individuo(1);

cout << i1->toString()<< " : " << i1->get_distancia() << endl;
cout << i2->toString()<< " : " << i2->get_distancia() << endl;
   vector <TIndividuo *> v = cruzamento->processa(i1, i2, populacao, getMaxGeracao(), 1);
cout << v[0]->toString()<< " : " << v[0]->get_distancia() << endl;
//cout << v[1]->toString()<<endl;
	
return;
*/

   populacao->ordena();
   time(&sysTime1);

   for (int i=1; i<=getMaxGeracao(); i++)
   {
      populacao = iteracao (populacao, i);
	  
      melhor = populacao->get_melhor();	  
	  if (melhor->get_distancia()<VP_melhor_dist)
	  {
         VP_melhor_dist = melhor->get_distancia();	
			VP_melhor_gera = i;
	  }
   }

   if (getPrintParcial())
   {
      cout << endl;
   }

   delete mutacao;
   delete cruzamento;

   time(&sysTime2);

   string strFimExec = "Temmpo de execução;";
   strFimExec += to_string(difftime(sysTime2, sysTime1));
   VP_ArqSaida->addLinha(strFimExec);
   if (getPrintParcial()) cout << strFimExec << endl;
	   
   strFimExec = "melhor;";
   strFimExec += to_string(melhor->get_distancia());
   VP_ArqSaida->addLinha(strFimExec);
   if (getPrintParcial()) cout << strFimExec << endl;
	   
   strFimExec = "Individuo;";
   strFimExec += melhor->toString();
   VP_ArqSaida->addLinha(strFimExec);
   if (getPrintParcial()) cout << strFimExec << endl;
	
   strFimExec = "Melhor Alcançado;";
   strFimExec += to_string(get_melhor_dst());
   strFimExec += "; Geração;";
   strFimExec += to_string(get_melhor_gera());
   VP_ArqSaida->addLinha(strFimExec);
   if (getPrintParcial()) cout << strFimExec << endl;

   strFimExec = "Qtde de Entradas Recursivas;";
   strFimExec += to_string(getEntradaRec());
   strFimExec += "; Qtde Execuções Recursivas;";
   strFimExec += to_string(getExecRec());
   strFimExec += "; Qtde Execuções das Combinações;";
   strFimExec += to_string(getCombinaRec());
   VP_ArqSaida->addLinha(strFimExec);

   VP_ArqSaida->addLinha("");
   
   delete populacao;
}

TPopulacao *TAlgGenetico::iteracao (TPopulacao *populacao, int geracao)
{
   TPopulacao *novaPop = new TPopulacao(populacao->get_tamanho(), populacao->getMapa(), populacao->getArqLog());
   elitismo(novaPop, populacao);
   geraDescendentes(novaPop, populacao, geracao);
   fazMutacao(novaPop);
   novaPop->ordena();

   delete populacao;

   time(&sysTime2);
   TIndividuo *melhor = novaPop->get_melhor();  
   TIndividuo *pior = novaPop->get_pior();

   if (getPrintParcial() && (get_profundidade()==0))
   {
      string lin = "Geração ";
      lin += to_string (geracao);
      lin += " / ";
      lin += to_string (getMaxGeracao());
      lin += ":";
   
      cout << lin << "(" << get_profundidade() << ") Melhor: " << melhor->get_distancia() << " - Pior: " << pior->get_distancia() << " - Média: " << novaPop->distanciaMedia() << " - Tempo: " << difftime(sysTime2, sysTime1) << endl;
      VP_ArqSaida->addTexto(lin);	  
   }

   VP_ArqSaida->addTexto(";Profundidade;");
   VP_ArqSaida->addTexto(to_string (get_profundidade()));
   VP_ArqSaida->addTexto(";Melhor;");
   VP_ArqSaida->addTexto(to_string (melhor->get_distancia()));
   VP_ArqSaida->addTexto(";Pior;");
   VP_ArqSaida->addTexto(to_string (pior->get_distancia()));
   VP_ArqSaida->addTexto(";Media;");
   VP_ArqSaida->addTexto(to_string(novaPop->distanciaMedia()));
   VP_ArqSaida->addTexto(";Tempo;");
   VP_ArqSaida->addLinha(to_string(difftime(sysTime2, sysTime1)));

   return novaPop;
}

void TAlgGenetico::elitismo(TPopulacao *novaPop, TPopulacao *populacao)
{
   TIndividuo *individuo;
   unsigned qtdeElitismo = populacao->get_tamanho() * getPercentElitismo() / 100;

   unsigned i = 0;
   double tamAnt = 0;

   while ((novaPop->get_qtdeIndividuo()<=qtdeElitismo)&&(i<populacao->get_qtdeIndividuo()))
   {
      individuo = populacao->get_individuo(i);
      if (individuo->get_distancia()!=tamAnt)
      {
         tamAnt = individuo->get_distancia();
         novaPop->add_individuo(individuo->clona());
      }
      i++;
   }
}

void TAlgGenetico::geraDescendentes(TPopulacao *novaPop, TPopulacao *populacao, int geracao)
{
   int i;
   int j;
   TIndividuo *ind1;
   TIndividuo *ind2;

   TSelecao *selecao = new TSelecao (VP_Mapa, VP_ArqSaida, VP_selecao);
	
   while (novaPop->get_qtdeIndividuo()<novaPop->get_tamanho())
   {
      i = selecao->processa(populacao);
      j = selecao->processa(populacao, i);

      ind1 = populacao->get_individuo(i);
      ind2 = populacao->get_individuo(j);

      if (ind1->get_distancia()==ind2->get_distancia())
      {
         ind1 = ind2->clona();
         novaPop->add_individuo(ind1);
      }
      else
		{
			novaPop->add_individuo(cruzamento->processa(ind1, ind2, populacao, getMaxGeracao(), geracao));
		}
         
   }
	
	delete selecao;
}

void TAlgGenetico::fazMutacao(TPopulacao *populacao)
{
   int ini=0;
   int qtdeMutacao          = populacao->get_tamanho() * getPercentMutacao() / 100;
	int qtdeMutacaoRecursiva = populacao->get_tamanho() * getPercentMutacaoRecursiva() / 100;

   switch (getSelIndMutacao())
   {
      case 0:  //O melhor nunca sofre mutação
      {
         ini = 1;
         break;
      }
      case 1:  //O melhor sempre sofre mutação
      {
         populacao->sub_dist_tot((populacao->get_individuo(0))->get_distancia());
         mutacao->processa(populacao->get_individuo(0));
         populacao->soma_dist_tot((populacao->get_individuo(0))->get_distancia());
         ini = 1;
         qtdeMutacao--;
         break;
      }
      default:  //A mutação do melhor indivíduo é aleatória
      {
         break;
      }
   }

   for (int i=0; i < qtdeMutacao; i++)
   {
      int j = TUtils::rnd (ini, populacao->get_tamanho()-1);

      populacao->sub_dist_tot((populacao->get_individuo(j))->get_distancia());

      mutacao->processa(populacao->get_individuo(j));

      populacao->soma_dist_tot((populacao->get_individuo(j))->get_distancia());
   }

   for (int i=0; i < qtdeMutacaoRecursiva; i++)
   {
      int j = TUtils::rnd (1, populacao->get_tamanho()-1);

      populacao->sub_dist_tot((populacao->get_individuo(j))->get_distancia());

      mutacaoAGRecursivo(populacao, j);

      populacao->soma_dist_tot((populacao->get_individuo(j))->get_distancia());
   }
	
}

void TAlgGenetico::mutacaoAGRecursivo(TPopulacao *populacao, int indice)
{
	VP_Entr_Rec++;

   //Tranosformará as sequências, do indivíduo, que coincidem com o melhor em um gene
   TTabConversao *tabConversao;

   int qtdeGenes;
   long esforco;
   long fatorial;

   TIndividuo *melhor = populacao->get_melhor();
   TIndividuo *manipulado = populacao->get_individuo(indice);

   //Verifico se o esforço
   //(poulação (cruzamento + elitismo))* Número de gerações * (população * %mutação)
   esforco = (populacao->get_tamanho() * getPercentMutacao()/100) * getMaxGeracao() * populacao->get_tamanho();

   //So o fatorial for <= 0 é porque houve overflow
   //-1, pois como é um ciclo, é possível fixar o 0
   fatorial = TUtils::fatorialStirling(manipulado->get_qtdeGenes()-1);

   //é maior do que todas as combinações (fatorial) da quantidade de gene
   //Se a combinação for maior, realiza o AG com os novos genes.
   //Obs.: o fatorial deve ser > 0, pois provavelmente irá estourar o long
   if ((fatorial>0)&&(esforco>=fatorial))
   {
   	VP_Comb_Rec++;
		
      //Quando o esforço é maior que as combinações
      //O melhor indivuo por combinação
		manipulado->melhorPossivel();

		if (manipulado->get_distancia()<melhor->get_distancia())
      {
         populacao->troca(0, indice);

			if(getPrintParcial()) cout << "Melhorou na combinação" << endl;
         VP_ArqSaida->addLinha ("Melhorou na combinação");
		}
      return;
   }

   do
   {
      // tentará reduzir a quantidade de gene.
      tabConversao = new TTabConversao(melhor, manipulado);
      qtdeGenes = tabConversao->get_quantidade();
	  
		//Obs.: qtdeGenes<=2, pois se for 1, indivíduo é igual ao melhor
      //      se for 2, não posso trocar as posições, já que sempre inicia
      //      pelo gene 0. Logo, posso sair sem fazer nada
      if (qtdeGenes<=2)
      {
         //Libera memória
         delete tabConversao;
         manipulado->inverte_sub_indice(TUtils::rnd(1, manipulado->get_qtdeGenes()-1), TUtils::rnd(1, manipulado->get_qtdeGenes()-1));
      }

   }while (qtdeGenes<=2);

   // Se não reduzir, não faz nada
	int reudcaoMinima = getPercentReducao()*0.01*manipulado->get_qtdeGenes();
   if (qtdeGenes >= reudcaoMinima)
   {
      //Libera memória
      delete tabConversao;
      return;
   }
   
	VP_Exec_Rec++;
   //Executando o AG recursivamente
   TAlgGenetico *ag = new TAlgGenetico(getMapa(), getArqLog(), get_profundidade()+1);
   ag->setMutacao(getMutacao());
   ag->setCruzamento(getCruzamento());
   ag->setTamPopulacao(getTamPopulacao());
   ag->setPrintParcial(getPrintParcial());
   ag->setMaxGeracao(getMaxGeracao());
   ag->setPercentElitismo(getPercentElitismo());
   ag->setPercentMutacao(getPercentMutacao());
   ag->setProfundidadeMaxima(getProfundidadeMaxima());
   ag->setSelecao(getSelecao());
   ag->setSelIndMutacao(getSelIndMutacao());
   ag->setPercentMutacaoRecursiva(getPercentMutacaoRecursiva());
   ag->setPercentReducao(getPercentReducao());
   ag->setTime(sysTime1);
   ag->exec(manipulado, tabConversao, melhor->get_extra());

   VP_Entr_Rec += ag->getEntradaRec();
	VP_Exec_Rec += ag->getExecRec();
	VP_Comb_Rec += ag->getCombinaRec();

   if (manipulado->get_distancia()<melhor->get_distancia())
   {
      populacao->troca(0, indice);
      string lin = "Melhorou profundidade ;";
	  
	  lin += to_string(get_profundidade()+1);
	  if (getPrintParcial()) cout << endl << lin << endl;
	  VP_ArqSaida->addLinha(lin);
   }
   else if (manipulado->get_distancia()==melhor->get_distancia())
   {
		//Realiza a troca de forma que mantenha uma redução
		int ind1 = TUtils::rnd(1, manipulado->get_qtdeGenes()-reudcaoMinima-1);
		int ind2 = ind1 + reudcaoMinima;
		
		//Garantindo que ind2 não será maior que a quantidade de genes
		if (ind2>=manipulado->get_qtdeGenes()) ind2 = manipulado->get_qtdeGenes()-1;
		
		if (ind2>ind1)
		   manipulado->inverte_sub_indice(TUtils::rnd(1, manipulado->get_qtdeGenes()-1), TUtils::rnd(1, manipulado->get_qtdeGenes()-1));
   }
   
   //Libera memória
   delete tabConversao;
}
