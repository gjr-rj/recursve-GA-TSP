#include "utils.hpp"
#include <iostream>

double TUtils::fatorialStirling (int n)
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

void TUtils::initRnd ()
{
   srand(time(NULL));
}

int TUtils::rnd(unsigned low, unsigned high)
{
   return rand() % (high-low+1) + low;
}

double TUtils::rndd(double low, double high)
{
	return rndd(low, high, 1000);
}

double TUtils::rndd(double low, double high, unsigned calibracao)
{
	//Gerando um número randomico entre 0 e VP_MaxrndD
   unsigned r = rand() % calibracao;
	
	//Após achar um randomico entre 0 e VP_MaxrndD
	//vou converter para a escala entre low e high
	return ((((double)r)/calibracao)*(high-low))+low;
}

/**
 *
 * Funcion que funciona como una moneda al aire, dependiendo de la
 * probabilidad que esta devuelve un 1 o un 0
 *
 **/
bool TUtils::flip(float prob)
{
   if (prob==0) return 0;
	
	double d = rnd(1, 10000);
	double p = prob * 10000;
  
  return (d <= p)?1:0;
}