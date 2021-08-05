#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' Interpolador/extrpolador linear
//'
//' Função para interpolar ou extrapolar linearmente em um conjunto de dados
//'
//' @param xData absissas dos dados
//' @param yData ordenadas dos dados
//' @param xPred absissas para interpolar/extrapolar
//' @param nExtrap número de observações utilizadas para extrapolação
//'
//' @return vetor numérico de dados interpolados/extrapolados
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector predictCpp2( NumericVector &xData, NumericVector &yData, NumericVector xPred, int nExtrap = 2)
{
  
  // Tamanho do vetor de previsao
  int sizePred = xPred.size();
  
  // Tamanho dos vetores de dados 
  int sizeData = xData.size();
  
  // Inicia vetor de saida
  NumericVector out(sizePred);
  
  // Checa se ult e maior que o numero de valores no dado. Se sim, o reduz para este numero
  if(nExtrap > sizeData)
  {
    nExtrap = sizeData;
  }
  
  // Inicializa variaveis da extrapolacao a direita (serao calculadas caso necessario)
  double betaDir = 0;
  double alfaDir = 0;
  
  // Testa se e preciso extrapolar a direita. Caso positivo, calcula os coeficientes
  if ( xPred[sizePred - 1] > xData[sizeData - 1])
  {
    // Extrai ultimos nExtrap elementos de xData e yData
    NumericVector xDataDir(nExtrap);
    
    for(int k = 0; k < nExtrap; ++k)
    {
      xDataDir[k] = xData[sizeData - nExtrap + k];
    }
    
    NumericVector yDataDir (nExtrap);
    
    for (int k = 0; k < nExtrap; ++k)
    {
      yDataDir[k] = yData[sizeData - nExtrap + k];
    }
    
    // Media desses ultimo nExtrap elementos
    double xMedDir = 0;
    double yMedDir = 0;

    for(int k = 0; k < nExtrap; ++k)
    {
      xMedDir += xDataDir[k] / nExtrap;
      yMedDir += yDataDir[k] / nExtrap;
    }
    
    // Calculo dos coeficientes da reta de extrapolacao
    double normDir = 0;
    for(int k = 0; k < nExtrap; ++k)
    {
      normDir += (xDataDir[k] - xMedDir) * (xDataDir[k] - xMedDir);
      betaDir += (xDataDir[k] - xMedDir) * (yDataDir[k] - yMedDir);
    }
    betaDir = betaDir / normDir;
    alfaDir = yMedDir - betaDir * xMedDir;
  }
  
  // Inicializa variaveis da extrapolacao a esquerda (serao calculadas caso necessario)
  double betaEsq = 0;
  double alfaEsq = 0;
  
  // Testa se e preciso extrapolar a esquerda. Caso positivo, calcula os coeficientes
  if ( xPred[0] < xData[0])
  {
    // Extrai primeiros nExtrap elementos de xData e yData
    NumericVector xDataEsq(nExtrap);
    
    for(int k = 0; k < nExtrap; ++k)
    {
      xDataEsq[k] = xData[k];
    }
    
    NumericVector yDataEsq (nExtrap);
    
    for (int k = 0; k < nExtrap; ++k)
    {
      yDataEsq[k] = yData[k];
    }
    
    // Media desses ultimo nExtrap elementos
    double xMedEsq = 0;
    double yMedEsq = 0;

    for(int k = 0; k < nExtrap; ++k)
    {
      xMedEsq += xDataEsq[k] / nExtrap;
      yMedEsq += yDataEsq[k] / nExtrap;
    }
    
    // Calculo dos coeficientes da reta de extrapolacao
    double normEsq = 0;
    for(int k = 0; k < nExtrap; ++k)
    {
      normEsq += (xDataEsq[k] - xMedEsq) * (xDataEsq[k] - xMedEsq);
      betaEsq += (xDataEsq[k] - xMedEsq) * (yDataEsq[k] - yMedEsq);
    }
    betaEsq = betaEsq / normEsq;
    alfaEsq = yMedEsq - betaEsq * xMedEsq;
  }
  
  // Loop para interpolar ou extrapolar o valor de cada x do vetor de previsao
  for(int i = 0; i < sizePred; ++i) {
    
    // Testa se o x em questao esta alem do dominio do ajuste e, se sim, extrapola
    if (xPred[i] < xData[0])
    {
      out[i] = alfaEsq + betaEsq * (xPred[i]);
    }
    else if ( xPred[i] > xData[sizeData - 1] )
    {
      out[i] = alfaDir + betaDir * (xPred[i]);
    }
    else
    {
      int j = 0;
      
      while ( xPred[i] > xData[j+1] ) j++;
      
      double xL = xData[j];
      
      double yL = yData[j];
      
      double xR = xData[j+1];
      
      double yR = yData[j+1];
      
      double dydx = ( yR - yL ) / ( xR - xL );
      
      out[i] = yL + dydx * ( xPred[i] - xL );
    }
    
  }
  
  return out;
}