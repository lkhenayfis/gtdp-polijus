######################## FUNCOES PARA CLASSIFICACAO E ANALISE DOS PATAMATRES #######################

#' Avaliação de patamares
#' 
#' Classifica o dado em patamares do nível de montante do aproveitamento a jusanta e os filtra
#' 
#' \\TODO: explicar como o processo de filtro funciona
#' 
#' @param dat objeto do tipo \code{datpoli} para classificar
#' @param tol vetor de tolerancias para o filro sucessivo. Ver Detalhes
#' @param plot.dir opcional. Se fornecido, diretorio onde salvar plots dos filtros aplicados a cada
#'     patamar -- pode ser fornecido \code{"auto"}; neste caso a raiz do arquivo de onde foram lidos
#'     os historicos e dados de extensao
#' 
#' @return objeto \code{datpoli} contendo item \code{infopats} preenchido e item \code{hist_est}
#'     filtrado de acordo
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' # leitura de uma planilha front
#' dat <- importadados()
#' 
#' # filtra vazoes estaveis
#' dat <- filtravazest(dat)
#' 
#' # classifica e filtra patamares
#' dat <- classfiltrapats(dat)
#' }
#' 
#' @export

classfiltrapats <- function(dat, tol = c(3, 2, 1.15), plot.dir) {

    # Classificacao
    hist_est <- dat$hist_est
    hist_est[, pat := round(nmont, 1)]
    hist_est[, pat := formatC(pat, digits = 1, format = "f", width = 5, flag = "0")]

    # limpa patamares individualmente
    hist_est <- split(hist_est, hist_est$pat)
    filtros <- lapply(hist_est, tratapats, tol = tol)
    hist_est <- rbindlist(hist_est)

    dat$hist_est <- hist_est
    dat$patinfo  <- filtros
    attr(dat, "classpats") <- TRUE

    if(!missing(plot.dir)) {
        if(plot.dir == "auto") {
            plot.dir <- attr(dat, "path")
            plot.dir <- strsplit(plot.dir, "\\\\")[[1]]
            plot.dir[length(plot.dir)] <- "Plots_Filtro"
            plot.dir <- do.call(file.path, as.list(plot.dir))
        }

        if(!dir.exists(plot.dir)) dir.create(plot.dir) else file.remove(list.files(plot.dir, full.names = TRUE))

        patamares <- unique(dat[[2]]$pat)

        for(pat in patamares) {

            out <- file.path(plot.dir, paste0("patamar_", pat, ".jpeg"))
            jpeg(out, , width = 2700, height = 1800, res = 300)
            plot(dat, paste0("pat_", pat))
            dev.off()
        }
    }

    return(dat)
}

# HELPERS ------------------------------------------------------------------------------------------

tratapats <- function(dpat, tol) {

    patamar <- dpat$pat[1]
    size <- nrow(dpat)

    niveisj <- dpat$njus

    #dpat[, elim := integer(size)]
    elim <- integer(size)

    # primeira eliminacao -- absurdos inferiores
    corta0 <- dpat[, min(njus)]
    if((size > 2) & (!all(niveisj == niveisj[1]))) {
        elim[niveisj <= quantile(niveisj, .03)] <- 1
    }

    tendencias <- lapply(seq(length(tol) + 1), function(x) NA)
    for(i in seq(tol)) {

        mod    <- tendencia(dpat[elim == 0])
        fitted <- if(class(mod) == "loess") mod$fitted else mod$fitted.values

        tendencias[[i]] <- mod

        # remove tendencia e calcula faixa de aceitacao
        dpat[elim == 0, njus := njus - fitted]
        iqr   <- dpat[elim == 0, diff(quantile(njus, c(.25, .75)))]
        corte <- dpat[elim == 0, range(quantile(njus, c(.25, .75)))] + tol[i] * c(-iqr, iqr)

        elim[(elim == 0) & ((dpat$njus < corte[1]) | (dpat$njus > corte[2]))] <-  i + 1

        dpat[elim %in% c(0, i + 1), njus := njus + fitted]
    }

    # ajusta tendencia final
    tendencias[[4]] <- tendencia(dpat[elim == 0])

    dpat[elim > 0, valido := FALSE]

    out <- list(tend = tendencias, filtro = elim)

    return(out)
}

tendencia <- function(d, span = 1) {
  if(nrow(d) > 20) {
    loess(njus ~ vazao, d, span = span, degree = 2)
  } else {
    lm(njus ~ vazao, d)
  }
}