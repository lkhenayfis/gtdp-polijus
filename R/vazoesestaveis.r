########################### FUNCOES PARA FILTRAGEM E TRATAMENTO DOS DADOS ##########################

#' Filtro de vazões estáveis
#' 
#' Função que realiza a filtragem de vazões estáveis no dado bruto
#' 
#' A identificação de vazões estáveis segue o processo descrito na 
#' \code{vignette("polijus", package = "polijus")}. De forma geral, são buscados \code{n} registros
#' horários consecutivos de vazão nos quais a máxima diferença percentual entre o primeiro e os 
#' \code{n - 1} restantes seja inferior ao valor \code{tol}.
#' 
#' A comparação é feita sempre em respeito ao primeiro registro e não de forma sequencial 
#' intencionalmente. Isto busca evitar a consideração incorreta de tomadas de carga lentas (ao longo
#' de múltiplas horas) como um registro estável. Por exemplo, se \code{n = 6} e \code{tol = .05}
#' 
#' \tabular{rr}{
#'     \eqn{Vazao_1} \tab 10.0\cr
#'     \eqn{Vazao_2} \tab 10.3\cr
#'     \eqn{Vazao_3} \tab 9.80\cr
#'     \eqn{Vazao_4} \tab 9.90\cr
#'     \eqn{Vazao_5} \tab 10.2\cr
#'     \eqn{Vazao_6} \tab 10.1
#' }
#' 
#' Todos os registros de 2 a 6 variam menos de 5% em relação ao primeiro, determinando assim uma
#' sequência de vazões estáveis. Caso estivéssmos numa condição de tomada de carga, com os mesmos 
#' parâmetros de execução
#' 
#' \tabular{rr}{
#'     \eqn{Vazao_1} \tab 10.0\cr
#'     \eqn{Vazao_2} \tab 10.3\cr
#'     \eqn{Vazao_3} \tab 10.7\cr
#'     \eqn{Vazao_4} \tab 11.1\cr
#'     \eqn{Vazao_5} \tab 11.6\cr
#'     \eqn{Vazao_6} \tab 12.0
#' }
#' 
#' Cada vazão está dentro da faixa de 5% daquela imediatamente anterior, porém, em função da subida
#' contínua, o último registro está a 20% de distância do primeiro. 
#' 
#' Uma vez determinadas todas as sequências de \code{n} vazões estáveis no dado, é obtido um 
#' histórico equivalente estável composto das médias de cada \code{n} vazões, níveis de jusante e 
#' níveis de montante do reservatório a jusante.
#' 
#' @param dat objeto \code{datpoli}. Ver \code{\link{new_datpoli}}
#' @param n inteiro indicando número de registros consecutivos buscados para sequência de vazões
#'     estáveis. Ver Detalhes
#' @param tol percentual no formato decimal (e.g. 0.1 para 10%) de variação máxima tolerada na
#'     sequência para determinação de estabilidade. Ver Detalhes
#' 
#' @return objeto \code{datpoli} com item \code{hist_est} preenchido com histórico equivalente
#'     estável. Ver Detalhes
#' 
#' @export

filtravazest <- function(dat, n = 6, tol = .05) {

    estavel <- NULL

    if(class(dat) != "datpoli") {
        stop("Parametro dat nao tem classe 'datpoli' -- Foi utilizada a funcao polijus::new_datpoli para contrui-lo?")
    }

    valido <- dat$hist$valido
    vazoes <- dat$hist$vazao
    tempo  <- as.numeric(dat$hist$datahora)
    size <- length(vazoes)

    # Inicializa janela de tamanho n e vetor de tag das sequencias estaveis no historico
    seqestavel <- rep(NA_real_, n)
    seqestavel[1] <- vazoes[1]
    tagestavel <- integer(size)

    for(i in 2:size) {

        if((diff(tempo[(i - 1):i]) != 3600) | !(valido[i])) {
            seqestavel <- rep(NA_real_, n)
            next
        }

        prox <- which(is.na(seqestavel))[1]
        seqestavel[prox] <- vazoes[i]

        consis <- abs(seqestavel - seqestavel[1]) / seqestavel[1]
        if(!all(consis[!is.na(consis)] < tol)) {

            # Caso haja mais de dois, vai removendo o mais a esquerda e testando de novo ate sobrar so um
            ok <- FALSE
            while(!ok) {

                seqestavel[1] <- NA_real_
                seqestavel <- shift(seqestavel, -1)

                if(sum(!is.na(seqestavel)) == 1) break

                consis <- abs(seqestavel - seqestavel[1]) / seqestavel[1]

                if(all(consis[!is.na(consis)] < tol)) ok <- TRUE
            }
        }

        # Uma vez que a janela estaja preenchida, marca as posicoes em tagestavel e reset
        if(all(!is.na(seqestavel))) {

            tagestavel[(i - n + 1):i] <- i
            seqestavel <- rep(NA_real_, n)
        }
    }

    # Toma medias do historico nas janelas estaveis
    hist_est <- copy(dat$hist)
    hist_est[, estavel := tagestavel]
    hist_est <- hist_est[estavel != 0, lapply(.SD, meanunique), by = estavel, .SDcols = 1:5]
    hist_est[, valido := as.logical(valido)]

    # Compoe saida
    dat$hist_est <- hist_est[, .SD, .SDcols = 2:6]
    attr(dat, "filtravazest") <- TRUE

    return(dat)
}

# HELPERS ------------------------------------------------------------------------------------------

shift <- function(vec, lag) {
    size <- length(vec)
    index <- seq(vec)
    if(sign(lag) == 1) {
        vec[c(tail(index, lag), head(index, size - lag))]
    } else {
        vec[c(tail(index, size - abs(lag)), head(index, abs(lag)))]
    }
}

# Esta funcao se faz necessaria pois em casos de media de n floats iguais usualmente tem um erro de
# arredondamento na ordem de 1e-15 ou menos. Por mais insignificante que pareca, isso pode levar a
# flutuacoes na classificacao de patamares.
# por exemplo, suponhamos a media de seis valores iguais a 445.450000000000000
# Em determinados casos a media deste vetor pode resultar em 445.449999999999999.
# O valor original, quando arrendodado para uma casa decimal vale 445.5, enquanto o segundo daria
# 445.4, gerando inconsistencia na classificacao.

meanunique <- function(vec, tol = 1e-5) {

    diferentes <- all(abs(vec - vec[1]) < tol)
    if(!diferentes) return(mean(vec)) else return(vec[1])
}
