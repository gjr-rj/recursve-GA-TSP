/*
*  arqlog.hpp
*
*  Módulo responsável pela geração do arquivo de log
*  Autor: Geraldo José Ferreira Chagas Junior - gjr.doc@gmail.com
*
*  PPGI - NCE - UFRJ
*  Data Criação: 06/05/2017
*  Datas de Modificações:
*
*/

#ifndef _ARQLOG_H
#define	_ARQLOG_H


#include <string>
#include <fstream>

using namespace std;

class TArqLog
{
   private:
      fstream fl;
      string VP_Path;

   public:
      TArqLog(string cabecalho, string path);
      ~TArqLog();
      void addLinha (string linha);
      void addTexto (string texto);
};

#endif
