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

#include "arqlog.hpp"

TArqLog::TArqLog(string cabecalho, string path)
{
   VP_Path = path;
   fl.open (VP_Path, std::fstream::out);
   addLinha (cabecalho);
}

TArqLog::~TArqLog()
{
   fl.close();
}
void TArqLog::addLinha (string linha)
{
   addTexto(linha);
   fl << endl;
}

void TArqLog::addTexto (string texto)
{
   fl << texto;
}
