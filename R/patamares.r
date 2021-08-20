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
#' # filtra vazoes estaveis
#' dat <- filtravazest(dummydata)
#' 
#' # classifica e filtra patamares
#' dat <- classfiltrapats(dat)
#' 
#' @export

classfiltrapats <- function(dat, tol = c(3, 2, 1.25), plot.dir) {

    nmont <- NULL

    if(!attr(dat, "filtravazest")) {
        stop("'dat' ainda nao passou pelo filtro de vazoes estaveis -- use polijus::filtravazest")
    }

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
    attr(dat, "classfiltrapats") <- TRUE

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

#' Extração de patamares
#' 
#' Função para extração de até quatro patamares para ajuste individual
#' 
#' A determinação de quais patamares sao ajustados através de curvas individuais se da visando 
#' encontrar um conjunto de curvas que, junto da curva base, estão o mais igualmente distribuídas 
#' pela dispersão de dados observados quanto possível. Atualmente o modelo DECOMP é capaz de lidar 
#' com até cinco famílias de polinômios -- contando a curva base, restam então no máximo quatro 
#' curvas associadas a patamares individuais do reservatório a jusante.
#' 
#' @param dat objeto do tipo \code{datpoli}
#' @param polibase objeto do tipo \code{polijusU} contendo ajuste da curva base
#' @param min_reg número mínimo de registros para eligibilidade a ajuste individual
#' @param quais vetor de strings indicando quais patamares extrair
#' 
#' @return lista de até quatro elementos contendo dados isolados de patamares
#' 
#' @export

extraipats <- function(dat, polibase, min_reg, quais) {

    dadospat <- split(dat$hist_est, dat$hist_est$pat)

    remanso   <- sapply(dat$patinfo, "[[", "remanso")
    dadospat  <- dadospat[remanso]

    if(!missing("quais")) {
        quais <- as.numeric(quais)
        quais <- formatC(quais, 1, 5, "f", "0")

        out <- dadospat[quais]

        if(any(sapply(out, function(d) is.null(dim(d))))) {
            qualna <- quais[which(sapply(out, function(d) is.null(dim(d))))]
            cli::cli_abort(c("patamares {qualna} nao podem ser utilizados",
                "x" = "Possiveis causas: patamar por inteiro nao tem remanso ou nao tem pontos"))
        }

        out <- lapply(out, function(d) d[valido == TRUE])

        return(out)
    }

    patamares <- names(dadospat)
    for(pat in patamares) {

        dpat <- dadospat[[pat]]
        vazconv <- dat$patinfo[[pat]]$vazconv_reg

        if(is.na(vazconv)) {

            mod <- tail(dat$patinfo[[pat]]$tend, 1)[[1]]

            vx  <- seq(dpat$vazao[1], dpat$vazao[nrow(dpat)], length.out = 500)
            vy1 <- predict(mod, data.frame(vazao = vx))
            vy2 <- predict(polibase, data.frame(vazao = vx))
            if(any(vy1 < vy2)) dpat[vazao > vx[which(vy1 < vy2)[1]], valido := FALSE]
        } else {

            dpat[vazao >= vazconv, valido := FALSE]
        }
    }

    numregpar <- sapply(dadospat, function(d) d[valido == TRUE, .N])
    patregmin <- patamares[numregpar > min_reg]
    patregminN <- as.numeric(patregmin)

    if(!(attr(dat, "namax") %in% as.numeric(patregmin))) {

        # // TODO:
        # implementar um metodo de selecao dos graus de liberdade e faixa para amostragem
        # do novo conjunto de dados

        #surf <- mgcv::gam(njus ~ s(vazao, k = 15, bs = "cs", m = 2) + s(nmont, k = 15, bs = "cs", m = 2),
        #    data = dat$hist_est[valido == TRUE])
    }

    numpoli <- ifelse(length(patregmin) < 4, length(patregmin), 4)

    minref   <- predict(polibase, newdata = data.frame(vazao = 0))
    rangeref <- seq(minref, tail(patregminN, 1), length.out = numpoli + 1)[-1]

    dists <- lapply(seq(numpoli), function(i) abs(patregminN - rangeref[i]))
    out   <- sapply(dists, function(v) patregmin[which.min(v)])

    out <- dadospat[out]

    out <- lapply(out, function(d) d[valido == TRUE])

    return(out)
}

# HELPERS ------------------------------------------------------------------------------------------

tratapats <- function(dpat, tol) {

    njus <- valido <- NULL

    size <- nrow(dpat)

    niveisj <- dpat$njus

    #dpat[, elim := integer(size)]
    elim <- integer(size)

    # primeira eliminacao -- absurdos inferiores
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
