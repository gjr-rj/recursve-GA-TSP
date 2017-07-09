/*
*  cruzamento.cpp
*
*  Módulo responsável pelas reotinas de cruzamento
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 19/05/2017
*  Datas de Modificações:
*
*/

#include "cruzamento.hpp"
#include <iostream>

//Apenas letura
TMapaGenes *TCruzamento::getMapa() { return VP_Mapa; }
TArqLog *TCruzamento::getArqLog()  { return VP_ArqSaida; }

TCruzamento::TCruzamento (TMapaGenes *mapa, TArqLog *arqSaida, int tipoCruzamento)
{
   VP_Mapa = mapa;
   VP_ArqSaida = arqSaida;
   VP_cruzamento = tipoCruzamento;
}

vector <TIndividuo *>TCruzamento::processa (TIndividuo *parceiro1, TIndividuo *parceiro2, TPopulacao *populacao, int maxGeracao, int geracao)
{
   switch (VP_cruzamento)
   {
      case 0:
      {
         return GSTX (parceiro1, parceiro2);
         break;
      }
      case 1:
      {
         return PMX(parceiro1, parceiro2);
         break;
      }
      case 2:
      {
			return OX1(parceiro1, parceiro2);
         break;
      }
      case 3:
      {
			return OX2(parceiro1, parceiro2);
         break;
      }
      case 4:
      {
			return MOX(parceiro1, parceiro2);
         break;
      }
      case 5:
      {
			return POS(parceiro1, parceiro2);
         break;
      }
      case 6:
      {
			return CX(parceiro1, parceiro2);
         break;
      }
      case 7:
      {			
			return DPX(parceiro1, parceiro2);
         break;
      }
      case 8:
      {
			return AP(parceiro1, parceiro2);
         break;
      }
      case 9:
      {			
			return MPX(parceiro1, parceiro2);
         break;
      }
      case 10:
      {
			return HX(parceiro1, parceiro2);
         break;
      }
      case 11:
      {
			return IO(parceiro1, populacao);
         break;
      }
      case 12:
      {
			return MIO(parceiro1, populacao, maxGeracao, geracao);
         break;
      }
      case 13:
      {
			return VR(parceiro1, parceiro2, populacao);
         break;
      }
      case 14:
      {
			return ER(parceiro1, parceiro2);
         break;
      }
      default:
      {
         break;
      }
   }

   vector <TIndividuo *> filhos;
   filhos.push_back(parceiro1->clona());
   filhos.push_back(parceiro2->clona());
   return filhos;
}

/**
 *
 * Edge Recombination Crossover (ER)
 *
 * Nguyen, Yoshihara, Yamamori & Yasunaga (2002). Greedy genetic algorithms
 * for symmetric and asymmetric TSPs. IPSJ Trans. Mathematical Modeling and
 * Its Applications, 43, 165-175.
 *
 * Sengoku & Yoshihara (1998). A fast TSP solver using GA on JAVA. In Third
 * International Symposium on Artificial Life, and Robotics (AROB III’98)
 * (pp. 283-288).
 *
 **/
vector <TIndividuo *>TCruzamento::GSTX(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
   TGene *gPar1;
   TGene *gPar2;

   vector<TGene *> temp;
   vector<bool> controle (parceiro1->get_qtdeGenes(), true);

   TIndividuo *filho = new TIndividuo(VP_Mapa, VP_ArqSaida);

   int pivo =  TUtils::rnd(1, parceiro1->get_qtdeGenes()-1);
   gPar1 = parceiro1->get_por_indice(pivo);
   gPar2 = parceiro2->get_por_id(gPar1->id);

   temp.push_back (gPar1);
   controle[gPar1->id] = false;

   bool dir = true;
   bool esq = true;
   int i=1;

   while (i<parceiro1->get_qtdeGenes())
   {
      if (esq)
      {
         if (gPar1->id==0) esq = false;
         else
         {
            gPar1 = parceiro1->ant(gPar1);
            if (controle[gPar1->id])
            {
               temp.insert(temp.begin(),gPar1);
               controle[gPar1->id] = false;
               i++;
            }
         }
      }

      if (dir)
      {
         if (gPar2->prox->id==0) dir = false;
         else
         {
            gPar2 = parceiro2->prox(gPar2);
            if (controle[gPar2->id])
            {
               temp.push_back (gPar2);
               controle[gPar2->id] = false;
               i++;
            }
         }
      }

      if ((!esq)&&(esq==dir))
      {
         gPar1 = parceiro1->ant(gPar1);
         if (controle[gPar1->id])
         {
            temp.push_back (gPar1);
            controle[gPar1->id] = false;
            i++;
         }

         gPar2 = parceiro2->prox(gPar2->id);
         if (controle[gPar2->id])
         {
            temp.push_back (gPar2);
            controle[gPar2->id] = false;
            i++;
         }
      }
   }

   filho->novo(temp);
   controle.clear();
   temp.clear();

   vector <TIndividuo *> filhos;
   filhos.push_back(filho);
   return filhos;
}

/**
 *
 * Partially-Mapped Crossover (PMX)
 *
 * Goldberg and Lingle (1985). Alleles, Loci and the TSP. In Grefenstette,
 * J. J. (ed.) Proceedings of the First International Conference on Genetic
 * Algorithms and Their Applications, 154–159. Hillsdale, New Jersey:
 * Lawrence Erlbaum.
 *
 **/
vector <TIndividuo *>TCruzamento::PMX(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
   int k;
   int uInicio, uFin;
   TGene *aux1;
   TGene *aux2;
   TIndividuo *filho1 = parceiro1->clona();
   TIndividuo *filho2 = parceiro2->clona();

   uInicio = TUtils::rnd(1, filho1->get_qtdeGenes()-2);
   uFin = TUtils::rnd(uInicio+1, filho1->get_qtdeGenes()-1);
  
   //O algorítimo está diferente do contreras, porém, ele foi
   //entendi e modificado de forma a ficar mais simples
   for(k=uInicio; k<=uFin; k++)
   {
      aux1 = filho1->get_por_indice(k);
      aux2 = filho2->get_por_indice(k);
      
      filho1->troca(aux1->id, aux2->id);
      filho2->troca(aux2->id, aux1->id);
   }
   vector <TIndividuo *> filhos;
   filhos.push_back(filho1);
   filhos.push_back(filho2);
   return filhos;
}

/**
 *
 * Order Crossover (OX1)
 *
 * Davis (1985). Applying Adaptive Algorithms to Epistatic Domains.
 * Proceedings of the International Joint Conference on Artiﬁcial Intelligence, 162–164
 *
 **/
vector <TIndividuo *>TCruzamento::OX1(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
   TIndividuo *filho2 = parceiro2->clona();

	filhos.push_back(filho1);
	filhos.push_back(filho2);

   //Para executar esse cruzamento a quantdade de genes tem 
	// que ser maior que 3.
	if (parceiro1->get_qtdeGenes()<=3) return filhos;
	
	int uInicio, uFin;
	uInicio = TUtils::rnd(1, filho1->get_qtdeGenes()-3);
	uFin = TUtils::rnd(uInicio+1, filho1->get_qtdeGenes()-2);

	TGene *it1;
	TGene *it2;
	int j = uFin;
	int fin_i = filho1->get_qtdeGenes();
	int fin_j = filho2->get_qtdeGenes();
	
	//Filho1
	for(int i = uFin+1; i < fin_i; i++)
	{
		while(j < fin_j)
		{
			j++;
         if(j >= parceiro2->get_qtdeGenes())
         {
				j = 1;
				fin_j = uFin;
         }
			
			it2 = parceiro2->get_por_indice(j);
			it1 = filho1->get_por_id(it2->id);
			if(it1->i<uInicio || it1->i>uFin) //Não está entre uInicio e uFin
         {
				filho1->troca_indice(i, it1->i);
				break;
         }
		}

      if(i >= filho1->get_qtdeGenes()-1)
      {
          i = 0;
          fin_i = uInicio;
      }
	}

	//Filho2
	j = uFin;
	fin_i = filho2->get_qtdeGenes();
	fin_j = filho1->get_qtdeGenes();
	for(int i = uFin+1; i < fin_i; i++)
	{
		while(j < fin_j)
		{
			j++;
         if(j >= parceiro1->get_qtdeGenes())
         {
				j = 1;
				fin_j = uFin;
         }
			
			it1 = parceiro1->get_por_indice(j);
			it2 = filho2->get_por_id(it1->id);
			if(it2->i<uInicio || it2->i>uFin) //Não está entre uInicio e uFin
         {
				filho2->troca_indice(i, it2->i);
				break;
         }
		}

      if(i >= filho2->get_qtdeGenes()-1)
      {
          i = 0;
          fin_i = uInicio;
      }
	}
	
	return filhos;
}

/**
 *
 * Order Based Crossover (OX2)
 *
 * Syswerda (1991). Schedule Optimization Using Genetic Algorithms. In Davis, L. (ed.)
 * Handbook of Genetic Algorithms, 332–349. New York: Van Nostrand Reinhold.
 *
 **/
vector <TIndividuo *>TCruzamento::OX2(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
   TIndividuo *filho2 = parceiro2->clona();

   // a quantidade de cidades a ser trocada
   unsigned qtdeTroca = parceiro1->get_qtdeGenes() * 0.4;

   //as posições selecionadas nos parceiros, já em ordem e sem repetição
   set<int> pPar1;
   set<int> pPar2;
			
   TGene *tmp;				
   int posicao;
	//Obs pPar1.size = pPar2.size. Se incluir em um, inclui no outro
	//Se não incluir em um, não incluirá no outro.
	while (pPar1.size() < qtdeTroca)
	{	
		//Seleção aleatória
		posicao = TUtils::rnd(1, parceiro1->get_qtdeGenes() - 1);
		pPar1.insert (posicao);
		
		//Identificando a cidade no parceiro 1
		tmp = parceiro1->get_por_indice(posicao);
		//Buscando a mesma cidade no Parceiro 2
		tmp = parceiro2->get_por_id(tmp->id);
		//Guardando a posição da cidade no Parceiro 2
		pPar2.insert(tmp->i);		
	}

   TGene *g1;
	TGene *g2;
	
	set<int>::iterator it1;
	set<int>::iterator it2;
   //Realizando as trocas
	for (it1=pPar1.begin(), it2 = pPar2.begin(); it1!=pPar1.end() && it2!=pPar2.end(); ++it1, ++it2)
	{
		//Troca no filho 1
		g1 = filho1->get_por_indice(*it1);
		g2 = parceiro2->get_por_indice(*it2);
		filho1->troca(g1->id, g2->id);
		
		//trocando no filho 2

	   g1 = filho2->get_por_indice(*it2);
		g2 = parceiro1->get_por_indice(*it1);
		filho2->troca(g1->id, g2->id);
	}

   pPar1.clear();
   pPar2.clear();
	
   filhos.push_back(filho1);
   filhos.push_back(filho2);
	
	return filhos;
}

/**
 *
 * Modified Order Based Crossover (MOX)
 *
 * New Operators of Genetic Algorithms for Traveling Salesman Problem
 *
 **/
vector <TIndividuo *>TCruzamento::MOX(TIndividuo *parceiro1, TIndividuo *parceiro2)
{

	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	
	filhos.push_back(filho1);

   //No mínimo os 2 primeiro genes serão mantidos ou 1 gene será trocado
	//Para isso, a quantidade de genes tem que ser no mínio 3
	if (parceiro1->get_qtdeGenes()<3) return filhos;

   int point = TUtils::rnd(2, parceiro1->get_qtdeGenes() - 1);
	
	TGene *tmp;
	for (int i=1, k=point; i < parceiro1->get_qtdeGenes(); i++)
	{
		//Pegando a cidade no parceiro 2
		tmp = parceiro2->get_por_indice(i);
		//Encontrando a mesma cidade no filho 1
      tmp = filho1->get_por_id(tmp->id);
		
		//Se a cidade estiver fora da parte imutável do filho
		//passa a fazer parte do filho na mesma sequência do parceiro 2
		if(tmp->i >= point)
		{
			filho1->troca_indice(k, tmp->i);
			k++;
		}
	}

   return filhos;
}

/**
 *
 * Position Based Crossover (POS)
 *
 * Syswerda (1991). Schedule Optimization Using Genetic Algorithms. In Davis, L. (ed.)
 * Handbook of Genetic Algorithms, 332–349. New York: Van Nostrand Reinhold.
 *
 **/
vector <TIndividuo *>TCruzamento::POS(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	
	filhos.push_back(filho1);

   // a quantidade de cidades a ser trocada
   unsigned qtdeTroca = parceiro1->get_qtdeGenes() * 0.4;

   //as posições no parceiro 1 a serem trocadas, já em ordem e sem repetição
   set<int> pPar1;
	//as posições dos genes que serão fixados 
	set<int> pPar2;
			
   int posicao;

	while (pPar1.size() < qtdeTroca)
	{	
		//Seleção aleatória
		posicao = TUtils::rnd(1, parceiro1->get_qtdeGenes() - 1);
		pPar1.insert (posicao);
	}

   TGene *g1;
	TGene *g2;
	
	set<int>::iterator it1;
	set<int>::iterator it2;

   //Realizando as trocas
	for (it1=pPar1.begin(); it1!=pPar1.end(); ++it1)
	{
		//Troca no filho 1 pelo referente na mesma posição
		//no parceiro 2
		g1 = filho1->get_por_indice(*it1);
		g2 = parceiro2->get_por_indice(*it1);
		filho1->troca(g1->id, g2->id);

      //Se g2 é fixado no filho guardo sua posição	Parceiro 1
		g1 = parceiro1->get_por_id(g2->id);
		pPar2.insert(g1->i);
	}
	
	//ajustando a ordem dos genes não fixados que devem
	//permanecer na mesma orde que parceiro 1
	it1 = pPar1.begin();
	it2 = pPar2.begin();
	
	for (int pos1=0, pos2=0; pos1 < parceiro1->get_qtdeGenes(); pos1++, pos2++)
	{
		//Aproveitando o fato de pPar1 e pPar2 estarem ordenados
		while (pos1==*it1) { pos1++; it1++; }
		while (pos2==*it2) { pos2++; it2++; }
		
		if (pos1 >= filho1->get_qtdeGenes() || pos2 >= parceiro1->get_qtdeGenes()) break;
		g1 = filho1->get_por_indice(pos1);
		g2 = parceiro1->get_por_indice(pos2);
		
		filho1->troca (g1->id, g2->id);
	}
	
   pPar1.clear();
   pPar2.clear();

	return filhos;
}

/**
 *
 * Cycle Crossover (CX)
 *
 * Oliver, Smith and Holland (1987). A Study of Permutation Crossover
 * Operators on the TSP. In Grefenstette, J. J. (ed.) Genetic Algorithms and
 * Their Applications: Proceedings of the Second International Conference,
 * 224–230. Hillsdale, New Jersey: Lawrence Erlbaum.
 *
 **/
vector <TIndividuo *>TCruzamento::CX(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	TIndividuo *filho2 = parceiro2->clona();
	
	filhos.push_back(filho1);
	filhos.push_back(filho2);

	//Selecionando o primeiro gene a ser trocado
	TGene *g1;
	TGene *g2;
	int ini;
	for (ini = 1; ini < filho1->get_qtdeGenes(); ini++)
	{
		g1 =  filho1->get_por_indice(ini);
		g2 =  parceiro2->get_por_indice(ini);
		
		if (g1->id != g2->id) break;
	}
	
   while (g1->id != g2->id)
	{
	   filho1->troca (g1->id, g2->id);
		filho2->troca (g1->id, g2->id);
		g2 = parceiro2->get_por_indice(g1->i);
	}
		
	return filhos;
}

/**
 *
 * Distance Preserving Crossover (DPX)
 *
 * Reisleben, Freisleben & Merz (1996). A genetic local search algorithm for
 * solving symmetric and asymmetric traveling salesman problems. Proceedings
 * of IEEE International Conference on Evolutionary Computation,
 * 1996. (pp. 616–621).
 *
 **/
vector <TIndividuo *>TCruzamento::DPX(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	filhos.push_back(filho1);

	vector<TGene *> fragmentoIni;
	vector<TGene *> fragmentoFim;
	TGene *g1;
	TGene *g2;
	
	bool primeiro = true;
	
	for (int i = 0; i < filho1->get_qtdeGenes(); i++)
	{
		g1 = filho1->get_por_indice(i);
		g2 = parceiro2->get_por_id(g1->id);
		
		if (primeiro) 
		{
			fragmentoIni.push_back(g1);
			primeiro = false;
		}
		
		if ((g1->prox->id != g2->prox->id)||(g1->prox->id==0))
      {
			fragmentoFim.push_back(g1);
         primeiro = true;
      }
	}

   int inverte;
	int proxFragmento;
	while(fragmentoIni.size() > 1)
	{
		double min = infinito, custo;
      int tmpInv;
	
      //Verificando quem tem menor custo
      for(unsigned i = 1; i < fragmentoIni.size(); i++)
		{
			g1 = fragmentoFim[0];
			g2 = fragmentoIni[i];
			double custo1 = VP_Mapa->get_distancia(g1->ori, g2->dest);
			g2 = fragmentoFim[i];
			double custo2 = VP_Mapa->get_distancia(g1->ori, g2->dest);

         if(custo1 <= custo2)
			{
				tmpInv = 0;
            custo = custo1;
			}
         else
         {
				tmpInv = 1;
				custo = custo2;
         }

         if(custo < min)
         {
				min = custo;
            inverte = tmpInv;
            proxFragmento = i;
			}
		}

		g1 = fragmentoFim[0];
		g2 = fragmentoIni[proxFragmento];
		
		int ini1 = g1->i + 1;
		int tam1 = 0;
		int ini2 = g2->i;
		
		g2 = fragmentoFim[proxFragmento];		
		int tam2 = g2->i - ini2 + 1;

      if(inverte) tam2 *= -1;
		
		filho1->troca_sub(ini1, tam1, ini2, tam2);
		
		if (inverte) fragmentoFim[0] = fragmentoIni[proxFragmento];
		else         fragmentoFim[0] = fragmentoFim[proxFragmento];

      fragmentoIni.erase(fragmentoIni.begin() + proxFragmento);
		fragmentoFim.erase(fragmentoFim.begin() + proxFragmento);
	}

   fragmentoIni.clear();
	fragmentoFim.clear();
	
	return filhos;
}

/**
 *
 * Alternating-position Crossover (AP)
 *
 * Larrañaga, Kuijpers, Poza & Murga (1997). Decomposing Bayesian networks:
 * triangulation of the moral graph with genetic algorithms. Statistics and
 * Computing 7, 1 (January 1997), 19-34.
 *
 **/
vector <TIndividuo *>TCruzamento::AP(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	TIndividuo *filho2 = parceiro2->clona();
	
	filhos.push_back(filho1);
	filhos.push_back(filho2);
	
	int iParc1 = 1;
	int iParc2 = 1;
   int iFilho1 = 1;
	int iFilho2 = 1;
	TGene *tmp;
	
	//Tratando o filho 1
   for(int i = 1; i < 2*parceiro1->get_qtdeGenes(); i++)
	{
		if(i%2)
		{
			if (iParc2>=parceiro2->get_qtdeGenes()) continue;
			//Buscando no parceito 2
			tmp = parceiro2->get_por_indice(iParc2);
			//Verificando no filho
			tmp = filho1->get_por_id(tmp->id);

			//Verificando se o Gene já foi incluido
			if (tmp->i == iFilho1) iFilho1++; //Ja está na posição
      	else if (tmp->i > iFilho1) //Ainda não está na posição correta
			{
				filho1->troca_indice(iFilho1, tmp->i);
				iFilho1++;
			}
			
			iParc2++;
		}
      else
      {
			if (iParc1>=parceiro1->get_qtdeGenes()) continue;
			//Buscando no parceito 2
			tmp = parceiro1->get_por_indice(iParc1);
			//Verificando no filho
			tmp = filho1->get_por_id(tmp->id);
			
			//Verificando se o Gene já foi incluido
			if (tmp->i == iFilho1) iFilho1++; //Ja está na posição
      	else if (tmp->i > iFilho1) //Ainda não está na posição correta
			{
				filho1->troca_indice(iFilho1, tmp->i);
				iFilho1++;
			}
			
			iParc1++;
      }
	}

	//Tratando o filho 2	
	iParc1 = 1;
	iParc2 = 1;
   for(int i = 1; i < 2*parceiro2->get_qtdeGenes(); i++)
	{
		if(i%2)
		{
			if (iParc1>=parceiro1->get_qtdeGenes()) continue;
			//Buscando no parceito 2
			tmp = parceiro1->get_por_indice(iParc1);
			//Verificando no filho
			tmp = filho2->get_por_id(tmp->id);
			
			//Verificando se o Gene já foi incluido
			if (tmp->i == iFilho2) iFilho2++; //Ja está na posição
      	else if (tmp->i > iFilho2) //Ainda não está na posição correta
			{
				filho2->troca_indice(iFilho2, tmp->i);
				iFilho2++;
			}
			
			iParc1++;
		}
      else
      {
			if (iParc2>=parceiro2->get_qtdeGenes()) continue;
			//Buscando no parceito 2
			tmp = parceiro2->get_por_indice(iParc2);
			//Verificando no filho
			tmp = filho2->get_por_id(tmp->id);
			
			//Verificando se o Gene já foi incluido
			if (tmp->i == iFilho2) iFilho2++; //Ja está na posição
      	else if (tmp->i > iFilho2) //Ainda não está na posição correta
			{
				filho2->troca_indice(iFilho2, tmp->i);
				iFilho2++;
			}
			
			iParc2++;
      }
	}

  return filhos;
}

/**
 *
 * Maximal preservative crossover (MPX)
 *
 * Mühlenbein, Gorges-Schleuter & Krämer (1988). Evolution algorithms in
 * combinatorial optimization. Parallel Computing, 7(1), 65–85.
 *
 * Mühlenbein (1991). Parallel genetic algorithms, population genetics and
 * combinatorial optimization. In J. D. Becker, I. Eisele, & F. W. Mündemann
 * (Eds.), Parallelism, Learning, Evolution (Vol. 565, pp. 398–406). Springer
 * Berlin Heidelberg.
 *
 **/
vector <TIndividuo *>TCruzamento::MPX(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	
	filhos.push_back(filho1);
	
	int begin, end;
	int Lmin, Lmax;

	Lmax = filho1->get_qtdeGenes()/2;

	if(filho1->get_qtdeGenes() <= 30) Lmin = TUtils::rnd(1, Lmax-1);
	else                              Lmin = 10;

	do 
	{
		begin = TUtils::rnd(1, filho1->get_qtdeGenes() - 2);
		end   = TUtils::rnd(begin + 1, filho1->get_qtdeGenes() - 1);
	}
	while((end - begin < Lmin) || (end - begin > Lmax));

   TGene *tmp;
	int index = 0;
   for (int i = 1; i < filho1->get_qtdeGenes(); i++)
	{
		//Primeiro, copia o bloco de parcero2 para o início do filho
	   if (begin <= end)	
		{
			//Quem é o indice no parceiro 2
			tmp = parceiro2->get_por_indice(begin);
			//Pegando o gene no filho
			tmp = filho1->get_por_id(tmp->id);
			filho1->troca_indice(i, tmp->i);
			begin++;
		}
		else  //pós preecher o início
		{
			do
			{
				index++;	
				
				//Não tem mais nada para ser posicionado
				if (index >= parceiro1->get_qtdeGenes()) break;
				
   			//Quem é o indice no parceiro 1
				tmp = parceiro1->get_por_indice(index);
				
				//Verficando se já está na posição correta, no filho
				tmp = filho1->get_por_id(tmp->id);
			}while (tmp->i < i);
			
			//Não tem mais nada para ser posicionado
			if (index >= parceiro1->get_qtdeGenes()) break;

			filho1->troca_indice(i, tmp->i);
		}		
	}
	
  return filhos;
}

/**
 *
 * Heuristic crossover (HX)
 *
 * Grefenstette, Gopal, Rosmaita & Van Gucht (1985). Genetic Algorithms for
 * the TSP. In Grefenstette, J. J. (ed.) Proceedings of the First
 * International Conference on Genetic Algorithms and Their Applications,
 * 160–165. Hillsdale, New Jersey: Lawrence Erlbaum.
 *
 **/
vector <TIndividuo *>TCruzamento::HX(TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	
	filhos.push_back(filho1);

	int i = TUtils::rnd(1, parceiro1->get_qtdeGenes() - 1);
	TGene *g = filho1->get_por_indice(i);
	TGene *found = parceiro1->get_por_id(g->id);

	TGene *gTmp;

	double custo;
	double min;
   TGene *proxG;

	//Obs.: por causa da forma que o indivíduo é tratado, tendo 0
	//sempre como a cidade inical, será necessário uma adequação.
	//O processo se fará até que encontre o zero. Então, o bloco
	//será movido para o final
	int ini = 1;
	int fim = filho1->get_qtdeGenes();
   filho1->troca_indice(ini, i);

	while(ini < fim-1)
   {
      min = infinito;

      //Verificando os vizinhos do gene em parente 1
		
		//se ainda não foi posicionado	
      gTmp = filho1->get_por_id(found->prox->id);	
		if ((gTmp->i > ini && gTmp->i < fim) || (gTmp->i==0 && fim == filho1->get_qtdeGenes()))
		{
		   custo = VP_Mapa->get_distancia(found->ori, found->prox->dest);
			if (custo < min)
			{
				min = custo;
				proxG = gTmp;
			}
		}
		
		//se ainda não foi posicionado		
      gTmp = filho1->get_por_id(found->ant->id);
		if ((gTmp->i > ini && gTmp->i < fim) || (gTmp->i==0 && fim == filho1->get_qtdeGenes()))
		{
		   custo = VP_Mapa->get_distancia(found->ori, found->ant->dest);
			if (custo < min)
			{
				min = custo;
				proxG = gTmp;
			}
		}

      //Verificando os vizinhos do gene em parente 2
		found = parceiro2->get_por_id(g->id);
			
		//se ainda não foi posicionado		
      gTmp = filho1->get_por_id(found->prox->id);
		if ((gTmp->i > ini && gTmp->i < fim) || (gTmp->i==0 && fim == filho1->get_qtdeGenes()))
		{
		   custo = VP_Mapa->get_distancia(found->ori, found->prox->dest);
			if (custo < min)
			{
				min = custo;
				proxG = gTmp;
			}
		}

		//se ainda não foi posicionado		
      gTmp = filho1->get_por_id(found->ant->id);	
		if ((gTmp->i > ini && gTmp->i < fim) || (gTmp->i==0 && fim == filho1->get_qtdeGenes()))
		{
		   custo = VP_Mapa->get_distancia(found->ori, found->ant->dest);     		
			if (custo < min)
			{
				min = custo;
				proxG = gTmp;
			}
		}

		//se nenhuma cidade esta disponível, sorteio
		if(min == infinito)
		{
         i = TUtils::rnd(ini+1, fim-1);
			proxG = filho1->get_por_indice(i);
		}
		
		//Tratando específicamente quando é o gene 0
		if (proxG->id == 0)
		{
			//Movo todo o bloco que já fo montado para o fim
			i = g->i;
			filho1->troca_sub(1, g->i, fim, 0);
			fim = fim - i;
			g = filho1->get_ini();
		}
		else
		{
			filho1->troca(g->prox->id, proxG->id);
   		g = g->prox;
		}
		
		found = parceiro1->get_por_id(g->id);
		ini = g->i;
	}

	return filhos;
}

/**
 *
 * Inver-over operator (IO)
 *
 * Tao & Michalewicz (1998). Inver-over operator for the TSP. In A. Eiben,
 * T. Bäck, M. Schoenauer, & H.-P. Schwefel (Eds.), Parallel Problem Solving
 * from Nature — PPSN V (Vol. 1498, pp. 803–812). Springer Berlin
 * Heidelberg.
 *
 **/
vector <TIndividuo *>TCruzamento::IO(TIndividuo *parceiro1, TPopulacao *populacao)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	
	float Prd = 0.02;
	
	int i = TUtils::rnd(1, parceiro1->get_qtdeGenes()-1);
	int j;
	TGene *c = filho1->get_por_indice(i);
	TGene *c_prim;
	
	while (true)
	{
		if(TUtils::flip(Prd))
		{
			do 
			{
            j = TUtils::rnd(1, filho1->get_qtdeGenes()-1);
            c_prim = filho1->get_por_indice(j);
			}
         while(c_prim->id == c->id);
		}
		else
		{
			j = TUtils::rnd(0, populacao->get_qtdeIndividuo()-1);
			TIndividuo *indSelect = populacao->get_individuo(j);
			c_prim = indSelect->get_por_id(c->id);
			c_prim = c_prim->prox;
         c_prim = filho1->get_por_id(c_prim->id);
			j = c_prim->i;
		}
		
		if ((c->prox->i == c_prim->i) || (c->ant->i == c_prim->i))
			break;

      
		if(i > j) filho1->inverte_sub(c->ant, c_prim);
      else      filho1->inverte_sub(c->prox, c_prim);
      
		c = c_prim;
		i = c->i;
	}
	
	if(parceiro1->get_distancia() < filho1->get_distancia())
	{
		delete filho1;
		filho1 = parceiro1->clona();		
	}
	
	filhos.push_back(filho1);
	return filhos;
}

/**
 *
 * Modified Inver-over operator (MIO)
 *
 * Wang, Sun, Li & Gao (2012). A Modified Inver-over Operator for the
 * Traveling Salesman Problem. In D.-S. Huang, Y. Gan, P. Gupta, &
 * M. M. Gromiha (Eds.), Advanced Intelligent Computing Theories and
 * Applications. With Aspects of Artificial Intelligence (Vol. 6839,
 * pp. 17–23). Springer Berlin Heidelberg.
 *
 **/
vector <TIndividuo *>TCruzamento::MIO (TIndividuo *parceiro1, TPopulacao *populacao, int maxGeracao, int geracao)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();

	float Prd = 0.02, PUCmax = 0.5, PUCmin = 0.2, Puc;
	int i = TUtils::rnd(1, parceiro1->get_qtdeGenes()-1);
	int j;
	unsigned qtdeMelhores = 5;
	TGene *c = filho1->get_por_indice(i);
	TGene *c_prim;
	
	while (true)
	{
		if(TUtils::flip(Prd))
		{
			vector<pair <int, double> > custos;
			for (int ii=1; ii<filho1->get_qtdeGenes(); ii++)
			{
				if (ii==c->i) continue;
				c_prim = filho1->get_por_indice(ii);
				VP_Mapa->get_distancia(c->ori, c_prim->dest);
				pair <int, double> custo (ii, VP_Mapa->get_distancia(c->ori, c_prim->dest));
				custos.push_back(custo);					
			}		
			sort (custos.begin(), custos.end(), compare_pair_second<std::less>());
							 
			j = TUtils::rnd(0, (custos.size()>qtdeMelhores)?qtdeMelhores-1:custos.size()-1);
			j = custos[j].first;
			c_prim = filho1->get_por_indice(j);
			
			custos.clear();
		}
		else
		{
			j = TUtils::rnd(0, populacao->get_qtdeIndividuo()-1);
			TIndividuo *indSelect = populacao->get_individuo(j);
			c_prim = indSelect->get_por_id(c->id);
			c_prim = c_prim->prox;
         c_prim = filho1->get_por_id(c_prim->id);
			j = c_prim->i;
		}
		
		if ((c->prox->i == c_prim->i) || (c->ant->i == c_prim->i))
			break;

      
		if(i > j)
		{
			int largo = filho1->get_qtdeGenes() - i + j, kkk, jjj;
			for(int kk = i + 1, r = 0, jj = j; r < (largo/2); r++, jj--, kk++)
			{
				if(jj == -1) jj = filho1->get_qtdeGenes() - 1;
				kkk = kk % filho1->get_qtdeGenes();
				jjj = jj % filho1->get_qtdeGenes();

				filho1->troca_indice (kkk, jjj);
			}
		}
      else      
			filho1->inverte_sub(c->prox, c_prim);

      Puc = PUCmax * exp((log(PUCmin / PUCmax)/maxGeracao)*geracao);
      if(TUtils::flip(Puc))
		{
			c = c_prim;
      }
		
		i = c->i;
	}
	
	if(parceiro1->get_distancia() < filho1->get_distancia())
	{
		delete filho1;
		filho1 = parceiro1->clona();		
	}
	
	filhos.push_back(filho1);
	return filhos;
}

/**
 *
 * Voting Recombination Crossover (VR)
 *
 * Mühlenbein (1989). Parallel Genetic Algorithms, Population Genetics and
 * Combinatorial Optimization. In Schaffer, J. (ed.) Proceedings on the Third
 * International Conference on Genetic Algorithms, 416–421. Los Altos, CA:
 * Morgan Kaufmann Publishers.
 *
 **/
vector <TIndividuo *>TCruzamento::VR (TIndividuo *parceiro1, TIndividuo *parceiro2, TPopulacao *populacao)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	filho1->embaralha();

	TSelecao *selecao = new TSelecao(VP_Mapa, VP_ArqSaida, TSelecao::tipoTorneioK);
	unsigned p3 = selecao->processa(populacao);
	unsigned p4 = selecao->processa(populacao);
	delete selecao;

	TIndividuo *parceiro3 = populacao->get_individuo (p3);
	TIndividuo *parceiro4 = populacao->get_individuo (p4);

   int qtdeIguais;
	int idGene;
	TGene *g;
	for (int i=1; i < filho1->get_qtdeGenes(); i++)
	{
		//Verifica se a posição i tem o mesmo gene em pelo menos 3
		// dos parceiros		
		g = parceiro1->get_por_indice(i);
		idGene = g->id;
		qtdeIguais = 1;

		g = parceiro2->get_por_indice(i);
		if (g->id==idGene) qtdeIguais++;
		
		g = parceiro3->get_por_indice(i);
		if (g->id==idGene) qtdeIguais++;
		
		g = parceiro4->get_por_indice(i);
		if (g->id==idGene) qtdeIguais++;
		
		if (qtdeIguais<3) 
		{
			g = parceiro2->get_por_indice(i);
			idGene = g->id;
			qtdeIguais = 1;

			g = parceiro3->get_por_indice(i);
			if (g->id==idGene) qtdeIguais++;
		
			g = parceiro4->get_por_indice(i);
			if (g->id==idGene) qtdeIguais++;			
		}
		
		if (qtdeIguais>=3) 
		{
			//Se encontrei um gene na mesma posição em pelomenos 3
			//parceiros, o filho, nesta posição, terá este gene
			g = filho1->get_por_indice(i);
			filho1->troca(g->id, idGene);
		}
	}

	filhos.push_back(filho1);
	return filhos;
}	

/**
 *
 * Edge Recombination Crossover (ER)
 *
 * Whitley, Starkweather & Fuquay (1989). Scheduling Problems and Travelling
 * Salesman: The Genetic Edge Recombination Operator. In Schaffer, J. (ed.)
 * Proceedings on the Third International Conference on Genetic Algorithms,
 * 133–140. Los Altos, CA: Morgan Kaufmann Publishers.
 *
 * Whitley, Starkweather & Shaner (1991). The Traveling Salesman and Sequence
 * Scheduling: Quality Solutions Using Genetic Edge Recombination. In Davis,
 * L. (ed.) Handbook of Genetic Algorithms, 350–372. New York: Van Nostrand
 * Reinhold.
 *
 **/
vector <TIndividuo *>TCruzamento::ER (TIndividuo *parceiro1, TIndividuo *parceiro2)
{
	vector <TIndividuo *> filhos;
	TIndividuo *filho1 = parceiro1->clona();
	
	TGene *g;
	TGene *gTemp;
	TGene *gVizinho;
	long custoTemp;
	
	int idMelhor;
	long custoMelhor;
	
	for (int i = 0; i < filho1->get_qtdeGenes()-2; i++ )
	{
		g = filho1->get_por_indice(i);
		
		//Servirá como a parte aleatória a medida que as posições	
		custoMelhor = VP_Mapa->get_distancia(g->ori, g->prox->dest);
      idMelhor = g->prox->id;
      
		//########################
		// Parceiro 1
		//########################
		//Pegando o mesmo gene no parcero 1
		gTemp = parceiro1->get_por_id(g->id);
		
		//Verificando se o próximo já está no filho
		//*****************************************
		gVizinho = filho1->get_por_id(gTemp->prox->id);
		if (gVizinho->i > i) //ainda não faz parte do filho
		{
			//Verifico o custo
			custoTemp = VP_Mapa->get_distancia(g->ori, gVizinho->dest);
			//Se tiver melhor custo, é um candidato a substituição
			if(custoTemp < custoMelhor)
			{
				custoMelhor = custoTemp;
				idMelhor = gVizinho->id;
			}
		}

		//Verificando se o anterior já está no filho
		//*******************************************
		gVizinho = filho1->get_por_id(gTemp->ant->id);
		if (gVizinho->i > i) //ainda não faz parte do filho
		{
			//Verifico o custo
			custoTemp = VP_Mapa->get_distancia(g->ori, gVizinho->dest);
			//Se tiver melhor custo, é um candidato a substituição
			if(custoTemp < custoMelhor)
			{
				custoMelhor = custoTemp;
				idMelhor = gVizinho->id;
			}
		}

		//########################
		// Parceiro 2
		//########################
		//Pegando o mesmo gene no parcero 2
		gTemp = parceiro2->get_por_id(g->id);
		
		//Verificando se o próximo já está no filho
		//*****************************************
		gVizinho = filho1->get_por_id(gTemp->prox->id);
		if (gVizinho->i > i) //ainda não faz parte do filho
		{
			//Verifico o custo
			custoTemp = VP_Mapa->get_distancia(g->ori, gVizinho->dest);
			//Se tiver melhor custo, é um candidato a substituição
			if(custoTemp < custoMelhor)
			{
				custoMelhor = custoTemp;
				idMelhor = gVizinho->id;
			}
		}

		//Verificando se o anterior já está no filho
		//*******************************************
		gVizinho = filho1->get_por_id(gTemp->ant->id);
		if (gVizinho->i > i) //ainda não faz parte do filho
		{
			//Verifico o custo
			custoTemp = VP_Mapa->get_distancia(g->ori, gVizinho->dest);
			//Se tiver melhor custo, é um candidato a substituição
			if(custoTemp < custoMelhor)
			{
				custoMelhor = custoTemp;
				idMelhor = gVizinho->id;
			}
		}

	   filho1->troca (g->prox->id, idMelhor);		
	}

	filhos.push_back(filho1);
	return filhos;
}	
