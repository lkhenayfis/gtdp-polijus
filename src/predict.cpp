#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
NumericVector predictCpp2( NumericVector &xData, NumericVector &yData, NumericVector xPred )
{
  
  // Tamanho do vetor de previsao
  int sizePred = xPred.size();
  
  // Tamanho dos vetores de dados 
  int sizeData = xData.size();
  
  // Inicia vetor de saida
  NumericVector out(sizePred);
  
  // Variavel de quantos pontos serao usados nas extrapolacoes
  int nExtrap = 2;
  
  // Checa se ult e maior que o numero de valores no dado. Se sim, o reduz para este numero
  if(nExtrap > sizeData)
  {
    nExtrap = sizeData;
  }
  
  // Inicializa variaveis da extrapolacao a direita (serao calculadas caso necessario)
  double Adir = 0;
  
  double Bdir = 0;
  
  // Testa se e preciso extrapolar a direita. Caso positivo, calcula os coeficientes
  if ( xPred[sizePred - 1] >= xData[sizeData - 1])
  {
    // Extrai ultimos cinco elementos de xData e yData
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
    
    // Cria vetor de xDataDir ao quadrado
    NumericVector xDataDir2 = xDataDir;
    
    for(int k = 0; k < nExtrap; ++k)
    {
      xDataDir2[k] = xDataDir2[k] * xDataDir2[k];
    }
    
    // Cria vetor de xDataDir * yDataDir
    NumericVector xyDataDir = xDataDir;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      xyDataDir[k] = xDataDir[k] * yDataDir[k];
    }
    
    // Realiza as somas dos vetores xDataDir, yDataDir, xyDataDir e xDataDir2
    double sumxDataDir = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumxDataDir += xDataDir[k];
    }
    
    double sumyDataDir = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumyDataDir += yDataDir[k];
    }
    
    double sumxyDataDir = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumxyDataDir += xDataDir[k] * yDataDir[k];
    }
    
    double sumxDataDir2 = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumxDataDir2 += xDataDir2[k];
    }
    
    // Resolve um minimos quadrados para extrapolar o dado
    Adir = (nExtrap * sumxyDataDir - sumxDataDir * sumyDataDir) / (nExtrap * sumxDataDir2 - sumxDataDir * sumxDataDir);
    
    Bdir = (sumyDataDir - Adir * sumxDataDir) / nExtrap;
  }
  
  // Inicializa variaveis da extrapolacao a esquerda (serao calculadas caso necessario)
  double Aesq = 0;
  
  double Besq = 0;
  
  // Testa se e preciso extrapolar a esquerda. Caso positivo, calcula os coeficientes
  if ( xPred[0] <= xData[0])
  {
    // Extrai primeiros cinco elementos de xData e yData
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
    
    // Cria vetor de xDataEsq ao quadrado
    NumericVector xDataEsq2 = xDataEsq;
    
    for(int k = 0; k < nExtrap; ++k)
    {
      xDataEsq2[k] = xDataEsq2[k] * xDataEsq2[k];
    }
    
    // Cria vetor de xDataEsq * yDataEsq
    NumericVector xyDataEsq = xDataEsq;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      xyDataEsq[k] = xDataEsq[k] * yDataEsq[k];
    }
    
    // Realiza as somas dos vetores xDataEsq, yDataEsq, xyDataEsq e xDataEsq2
    double sumxDataEsq = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumxDataEsq += xDataEsq[k];
    }
    
    double sumyDataEsq = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumyDataEsq += yDataEsq[k];
    }
    
    double sumxyDataEsq = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumxyDataEsq += xDataEsq[k] * yDataEsq[k];
    }
    
    double sumxDataEsq2 = 0;
    
    for (int k = 0; k < nExtrap; ++k)
    {
      sumxDataEsq2 += xDataEsq2[k];
    }
    
    // Resolve um minimos quadrados para extrapolar o dado
    Aesq = (nExtrap * sumxyDataEsq - sumxDataEsq * sumyDataEsq) / (nExtrap * sumxDataEsq2 - sumxDataEsq * sumxDataEsq);
    
    Besq = (sumyDataEsq - Aesq * sumxDataEsq) / nExtrap;
  }
  
  // Loop para interpolar ou extrapolar o valor de cada x do vetor de previsao
  for(int i = 0; i < sizePred; ++i) {
    
    // Testa se o x em questao esta alem do dominio do ajuste e, se sim, extrapola
    if (xPred[i] <= xData[0])
    {
      out[i] = Besq + Aesq * (xPred[i]);
    }
    else if ( xPred[i] >= xData[sizeData - 1] )
    {
      out[i] = Bdir + Adir * (xPred[i]);
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