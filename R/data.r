##################################### DADOS PARA EXEMPLIFICACAO ####################################

#' Exemplo de objeto \code{datpoli}
#' 
#' Exemplo de instância da classe \code{datpoli} obtida a partir da leitura da planilha front
#' dummydata.xlsx, também inclusa no pacote
#' 
#' @format
#' Lista de três elementos, sendo o primeiro o historico horario
#' \describe{
#'     \item{datahora}{POSIXt indicando momento da observacao}
#'     \item{njus}{valor do nível de jusante na data e hora especificada}
#'     \item{vazao}{vazão defluente total (turbinamento + vertimento) da usina}
#'     \item{nmont}{nível de montante do reservatório a jusante}
#'     \item{valido}{booleano indicando a usabilidade do registro horario}
#' }
#' 
#' O segundo elemento contem uma lista de dados de extensão. Neste caso há apenas uma fonte

"dummydata"