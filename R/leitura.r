############################# LEITURA E FORMATACAO DOS DADOS DE ENTRADA ############################

#' Leitura da planilha padrão
#' 
#' Função para leitura e preparo dos dados a partir da planilha front
#' 
#' Esta função é uma facilitadora para leitura dos dados de vazão e nível de jusante a partir da 
#' planilha padrão para entrada de dados. Apenas uma facilitadora pois é possível alcançar o mesmo
#' objetivo lendo os dados necessários a partir de qualquer outra fonte e local -- ver 
#' \code{\link{new_datpoli}} para mais informações.
#' 
#' @param path o caminho da planilha de entrada. Caso não seja fornecido, uma janela do explorador
#'     de arquivos será aberta para busca e seleção
#' 
#' @return objeto \code{datpoli} contendo dados historicos e funções para extensão dos mesmos
#' 
#' @examples 
#' 
#' # Fornecendo um caminho diretamente a funcao
#' dummypath <- system.file("extdata", "dummydata.xlsx", package = "polijus")
#' dat <- importadados(path = dummypath)
#' 
#' \dontrun{
#' # Chamando a funcao sem parametro path para abrir o explorador de arquivos
#' dat <- importadados()
#' }
#' 
#' @seealso \code{\link{new_datpoli}} para construtor do objeto de saída
#' 
#' @import data.table
#' @importFrom readxl read_xlsx
#' 
#' @export

importadados <- function(path = NULL) {

    if(is.null(path)) path <- file.choose()

    # Le cadastro e pega codigo da usina
    cad <- read_xlsx(path, sheet = 1)
    cod <- as.numeric(cad[1, 2])

    # Leitura do historico
    hist <- setDT(read_xlsx(path, sheet = 2))
    colnames(hist) <- c("data", "hora", "njus", "vazao", "nmont", "valido")

    # Leitura dos dados de extensao
    ext <- setDT(read_xlsx(path, sheet = "Extensao", col_types = "numeric", skip = 1))
    ext <- lapply(seq(ncol(ext))[c(TRUE, FALSE)], function(i) {
        cols <- c(i, i + 1)
        out <- ext[, ..cols]
        colnames(out) <- c("vazao", "njus")
        out
    })

    new_datpoli(path, cod, hist, ext)
}

#' Inicializador de objeto datpoli
#' 
#' Função para criação de novas instâncias de objetos da classe datpoli
#' 
#' Tipicamente esta função só é chamada ao final de \code{\link{importadados}} para consolidar as
#' informações lidas na planilha padrão em um único objeto de saída. Porém, como mencionado na 
#' documentação da mesma, tais objetos não precisam ser criados necessariamente fazendo uso da 
#' planilha padrão, sendo possível lê-los de outros locais e combiná-los com esta função. Para que 
#' isto ocorra sem erros, é necessário garantir algumas características dos dados.
#' 
#' Quanto aos dados históricos, estes devem ser um objeto em formato tabelado (data.frame, 
#' data.table, tibble etc) contendo colunas
#' \describe{
#'     \item{data}{data da observação em formato \code{Date}}
#'     \item{hora}{hora da observação como inteiro de 1 a 24}
#'     \item{njus}{valor do nível de jusante na data e hora especificada}
#'     \item{vazao}{vazão defluente total (turbinamento + vertimento) da usina}
#'     \item{nmont}{nível de montante do reservatório a jusante}
#'     \item{valido}{booleano indicando a usabilidade do registro horario}
#' }
#' 
#' Em segundo lugar, o parâmetro \code{ext} deve ser fornecido como uma lista de objetos tabulados
#' contendo as colunas
#' \describe{
#'     \item{vazao}{vazão defluente total}
#'     \item{njus}{nível de jusante}
#' }
#' 
#' Em ambos os casos as colunas podem estar fora de ordem, contanto que todas existam com os nomes
#' corretos.
#' 
#' @param cod escalar inteiro indicando o código da usina
#' @param hist dado histórico horário da usina. Ver Detalhes
#' @param ext lista de funções para extensão do dado histórico. Ver Detalhes
#' 
#' @return objeto \code{datpoli} contendo dados historicos e funções para extensão dos mesmos
#' 
#' @seealso \code{\link{importadados}} para leitura das informações a partir da planilha padrão.
#' 
#' @examples
#' 
#' # Dado historico aleatorio
#' d_hist <- data.frame(vazao = rnorm(240), nmont = rep(seq(10.1, 12, by = .1), each = 12),
#'     njus = rnorm(240), hora = 1:24, data = as.Date("2020-01-01") + rep(1:10, each = 24),
#'     valido = TRUE)
#' 
#' # Extensao
#' l_ext <- list(data.frame(vazao = 1:6, njus = log(1:6)), data.frame(vazao = 0:7, njus = .8 * (0:7) + 1))
#' 
#' # Usando codigo 6 (Furnas)
#' new_datpoli(cod = 6, hist = d_hist, ext = l_ext)
#' 
#' @export

new_datpoli <- function(path, cod, hist, ext) {

    hist <- tratahist(hist)
    ext  <- trataext(ext)

    vazef <- HIDR[HIDR$USI == cod, "VAZEF"]

    coef <- unlist(HIDR[HIDR$USI == cod, -1])[1:5]
    fun  <- function(x) coef[1] + coef[2] * x + coef[3] * x^2 + coef[4] * x^3 + coef[5] * x^4
    ext  <- c(list(CAD = fun), ext)

    new <- list(hist = hist, hist_est = NULL, ext = ext, patinfo = NULL)
    class(new) <- "datpoli"

    attr(new, "path")  <- path
    attr(new, "cod")   <- cod
    attr(new, "vazef") <- vazef
    attr(new, "estavel")   <- FALSE
    attr(new, "classpats") <- FALSE

    return(new)
}

#' @export

print.datpoli <- function(dat) print(dat$hist)

# HELPERS ------------------------------------------------------------------------------------------

#' @import data.table

tratahist <- function(hist) {

    if(!is.data.table(hist)) hist <- as.data.table(hist)

    colunas <- c("data", "hora", "njus", "vazao", "nmont", "valido")
    if(!all(colunas %in% colnames(hist))) {
        stop("Faltam as colunas\n", colnames(hist)[!(colunas %in% colnames(hist))])
    } else {
        setcolorder(hist, colunas)
    }

    hist[, datahora := paste(data, paste0((hora - 1), ":00:00"))]
    hist[, datahora := fasttime::fastPOSIXct(datahora, tz = "GMT")]
    setorder(hist, "datahora")
    if(!all(diff(hist$datahora) == 1)) stop("Ha datas faltantes no historico")

    colunas <- c("datahora", "njus", "vazao", "nmont", "valido")
    hist <- hist[, ..colunas]

    hist[vazao < 0, valido := FALSE]
    hist[njus < 0, valido := FALSE]
    hist[nmont < 0, valido := FALSE]

    hist[!complete.cases(hist), valido := FALSE]

    return(hist)
}

#' @import data.table

trataext <- function(ext) {

    if(class(ext) != "list") ext <- list(ext)
    for(i in seq(ext)) {
        ext[[i]] <- as.data.table(ext[[i]])
        if(!all(c("vazao", "njus") %in% colnames(ext[[i]]))) stop("Extensao #", i, " nao contem colunas necessarias")
    }
    if(length(ext) > 0) names(ext) <- paste0("EXT", seq(ext))

    ext <- ext[!sapply(ext, function(d) all(is.na(d)))]

    ext <- lapply(ext, function(d) {
        d   <- setorder(d, "vazao")
        mod <- loess(njus ~ vazao, d)
        function(x) predictCpp2(d$vazao, mod$fitted, x)
    })

    return(ext)
}