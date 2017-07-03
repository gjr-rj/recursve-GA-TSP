# recursve-GA-TSP
Recursive Genetic Algorithm for TSP resolver (in C++)

# compilação

./make

# linha de execução

./tsp <arquivo de instância XML> <arquivo de configuração> <inicio do nome do arquivo de saida a ser criado>

ex.:
./tsp files/gr96.xml files/tspconfig.xml resultados/r-gr96

onde:
* **./tsp** - Algoritmo genético recursivo;
* **files/gr96.xml** - arquivo de instancia XML com 96 cidades;
* **files/tspconfig.xml** - arquivo de confguração XM;
* **resultados/r-gr96** -  inicio do nome do arquivo de saida a ser criado.

exemplo de arquivos a serem criados são:
   - r-gr96-0-2017-07-01-10-23-20.txt;
   - r-gr96-1-2017-07-01-10-23-31.txt;
   - r-gr96-2-2017-07-01-10-23-50.txt,
                   
    onde :
      r-gr96 - parte do nome informado na execição;
      1, 2, 3 - resultados da primeira, segunda e terceira execução;
      2017-07-01 - data da execução;
      10-23-20, 10-23-31, 10-23-50 - hora, minuto e segundo de cada excução.
                       
# parametros de configuração (arquivo XML)

                       
 * **integer** - tamanho em digitos do inteiro;
 * **printParcial**
      - 1 imprime resultados parciais na tela,
      - 0 desativa  impressão;               
* **tamanhoPopulacao** - tamanho da população;
* **numGeracoes** - número de gerações;
* **percentMutacao** - perentual de mutação;
* **mutacao**
  - 0 Exchange mutation (EM)
  - 1 2-opt,
  - 2 2-opt (2),
  - 3 Fast 3-opt,
  - 4 3-opt,
  - 5 Scramble mutation (SM),
  - 6 Simple Inversion mutation (SIM),
  - 7 Displacement mutation (DM),
  - 8 Inversion mutation (IVM),
  - 9 Inversion mutation (ISM),
  - 10 Greedy-swap mutation (GSM),
  - 11 Heuristic mutation (HM),
  - 12 Greedy Sub Tour Mutation (GSTM),
  - 13 DBM,
  - 14 Self-adaptive Hybrid Mutation Operator (SHMO),
  - 15 Double Bridge move (DBM),
  - 16 TIPO3,
  - 17 TIPO4,
  - 18 Neighbor-Join (NJ);       
* **selIndMutacao** 
  - 0 o melhor indivíduo nunca sofrerá mutação,
  - 1 o melhor indivíduo sempre sofrerá mutação,
  - 2 a escolha do indivíduo para mutação é completamente aleatória;                
* **cruzamento**
  - 0 GSTX,
  - 1 Partially-Mapped Crossover (PMX),
  - 2 Order Crossover (OX1),
  - 3 Order Based Crossover (OX2),
  - 4 Modified Order Based Crossover (MOX),
  - 5 Position Based Crossover (POS),
  - 6 Cycle Crossover (CX),
  - 7 Distance Preserving Crossover (DPX),
  - 8 Alternating-position Crossover (AP),
  - 9 Maximal preservative crossover (MPX),
  - 10 Heuristic crossover (HX),
  - 11 Inver-over operator (IO),
  - 12 Modified Inver-over operator (MIO),
  - 13 Voting Recombination Crossover (VR),
  - 14 Edge Recombination Crossover (ER);
* **numExec** - quantidade de execuções (repetições)
* **roleta**
  - 0 desativa roleta,
  - 1 ativa roleta;
* **percentElitismo** - percentual de elitismo;
* **percentMutacaoRecursiva** - percentual de mutação recursiva;
* **percentReducao** - percentual de redução, na mutação recursiva;
* **profundidadeMaxima** - profundidade máxima na recursivdade;

# exemplo de arquivo de configuração
```
<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
<configuracaoTSP>

  <name>tspconfig</name>

  <source>TSPLIB</source>

  <description>Configuração da programa de Algoritmo Genético</description>

  <integer>15</integer>

  <printParcial>1</printParcial>
  <tamanhoPopulacao>300</tamanhoPopulacao>
  <numGeracoes>200</numGeracoes>
  <percentManipulacao>30</percentManipulacao>
  <percentMutacao>1</percentMutacao>
  <percentReducao>78</percentReducao>
  <mutacao>6</mutacao>
  <selIndMutacao>0</selIndMutacao>
  <cruzamento>12</cruzamento>
  <numExec>1</numExec>
  <roleta>0</roleta>
  <percentElitismo>70</percentElitismo>
  <percentMutacaoRecursiva>20</percentMutacaoRecursiva>
  <profundidadeMaxima>10</profundidadeMaxima>
</configuracaoTSP>```
