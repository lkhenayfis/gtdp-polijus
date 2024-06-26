############################ FUNCOES PARA AVALIACAO DO EFEITO DE REMANSO ###########################

#' Avaliacao do efeito de remanso
#' 
#' Identifica quais patamares sofrem influencia parcial ou completa do reservatorio a jusante
#' 
#' //TODO explicar melhor o processo de avaliacao de remanso
#' 
#' @param dat objeto do tipo \code{datpoli} ja classificado no qual avaliar remanso
#' @param tol tolerancia minima para consideracao de convergencia entre dois patamares. Ver Detalhes
#' @param step separacao entre patamares para testar covergencia. Ver Detalhes
#' 
#' @return objeto \code{datpoli} contendo item \code{infopats} preenchido com vazoes de convergencia
#' 
#' @export 

evalremanso <- function(dat, tol = .05, step = 2) {

    vazao <- njus <- NULL

    patamares <- sort(unique(dat$hist_est$pat))

    # Heuristica para consideracao de convergencia imediata dos patamares
    vazmin <- quantile(dat$hist_est[, vazao], .05) + diff(quantile(dat$hist_est[, vazao], c(.05, 1))) / 25

    # ANALISE DE CONVERGENCIA

    converg <- lapply(seq(patamares)[-c(1:step)], function(i) {

        pat1 <- patamares[i - step]
        pat2 <- patamares[i]

        checaconv(pat1, pat2, tol, vazmin, dat)
    })
    converg <- c(lapply(seq(step), function(x) c(vazconv = NA, convimed = NA)), converg)

    dat$patinfo <- mapply(dat$patinfo, converg, SIMPLIFY = FALSE, FUN = function(info, conv) {
        info$vazconv  <- unname(conv[1])
        info$convimed <- unname(conv[2])
        info
    })

    # ANALISE DE PRESENCA DE REMANSO

    # Para identificacao da ausencia de influencia de remanso, um dos testes e verificar se
    # o maior patamar esta significativamente abaixo do menor nivel de jusante da usina em questao
    # Este "significativamente" varia caso a caso, de modo que e necessaria uma heuristica geral
    # para realizacao deste teste

    iqr    <- diff(range(quantile(dat$hist_est[vazao < vazmin, njus], c(0.25, 0.75))))
    limjus <- quantile(dat$hist_est[vazao < vazmin, njus], 0.25) - 1.2 * iqr
    limjus <- min(dat$hist_est[(vazao < vazmin) & (njus > limjus), njus])

    semremanso <- achasemremanso(dat, vazmin, limjus)

    dat$patinfo <- mapply(dat$patinfo, semremanso, SIMPLIFY = FALSE, FUN = function(info, srms) {
        info$remanso <- !srms
        info
    })

    # REGULARIZACAO DAS VAZOES DE CONVERGENCIA

    vazconvreg <- tratavazconv(dat)
    dat$patinfo <- mapply(dat$patinfo, vazconvreg, SIMPLIFY = FALSE, FUN = function(info, vcr) {
        info$vazconv_reg <- vcr
        info
    })

    attr(dat, "step_converg") <- step
    attr(dat, "evalremanso") <- TRUE

    return(dat)
}

# HELPERS ------------------------------------------------------------------------------------------

checaconv <- function(pat1, pat2, tol, vazmin, dat) {

    vazao <- pat <- valido <- NULL

    d1 <- dat$hist_est[(pat == pat1) & (valido == TRUE)]
    d2 <- dat$hist_est[(pat == pat2) & (valido == TRUE)]

    m1 <- tail(dat$patinfo[[pat1]]$tend, 1)[[1]]
    m2 <- tail(dat$patinfo[[pat2]]$tend, 1)[[1]]

    if((max(d1$vazao) < min(d2$vazao)) | (max(d2$vazao) < min(d1$vazao))) {
        predx <- NA
    } else {
        v_range <- c(max(min(d1$vazao), min(d2$vazao)), min(max(d1$vazao), max(d2$vazao)))
        predx <- c(d1[(vazao > v_range[1]) & (vazao < v_range[2]), vazao],
                    d2[(vazao > v_range[1]) & (vazao < v_range[2]), vazao])
        predx <- predx[order(predx)]
    }
    if(length(predx) == 0) predx <- NA

    xdata1 <- if(class(m1) == "loess") m1$x[, 1] else m1$model[, 2]
    xdata2 <- if(class(m2) == "loess") m2$x[, 1] else m2$model[, 2]
    ydata1 <- if(class(m1) == "loess") m1$fitted else m1$fitted.values
    ydata2 <- if(class(m2) == "loess") m2$fitted else m2$fitted.values

    order1 <- order(xdata1)
    order2 <- order(xdata2)

    predy1 <- predictCpp2(xdata1[order1], ydata1[order1], predx)
    predy2 <- predictCpp2(xdata2[order2], ydata2[order2], predx)

    dif    <- predy2 - predy1
    diftol <- !is.na(dif) & (abs(dif) < tol)

    out <- c(vazconv = NA_real_, convimed = NA_real_)

    if(any(diftol)) {

        # Caso positivo, reduz o vetor de atendimento e vazoes apenas para do primeiro T em diante
        predx   <- predx[which(diftol)[1]:length(diftol)]
        diftol  <- diftol[which(diftol)[1]:length(diftol)]

        fim <- FALSE
        while(!fim) {

            # Calcula o indice de atendimento a tolerancia a partir da primeira vazao T
            indatend <- mean(diftol)

            # Checa se atende ao criterio de persistencia
            if(indatend > 0.9) {

                # Caso positivo, pega a vazao na qual a convergencia ocorre
                out[1] <- head(predx[diftol], 1)

                # Testa se a convergencia foi imediata
                if(out[1] <= vazmin) out[2] <- TRUE

                fim <- TRUE

            } else {

                # Caso negativo, reduz o vetor para o primeiro T apos o trecho F
                predx  <- predx[which(!diftol)[1]:length(diftol)]
                diftol <- diftol[which(!diftol)[1]:length(diftol)]

                tryCatch({
                    predx  <- predx[which(diftol)[1]:length(diftol)]
                    diftol <- diftol[which(diftol)[1]:length(diftol)]
                }, error = function(e) fim <<- TRUE)
            }
        }
    }

    return(out)
}

achasemremanso <- function(dat, vazmin, limjus) {

    vazao <- pat <- valido <- NULL

    pats <- sort(unique(dat$hist_est$pat))

    v_convimed <- sapply(dat$patinfo, "[[", "convimed") == 1

    # Outra heuristica para esse corte brusco e o - 3, significando um nivel maximo do reserv a
    # jusante muito menor que o minimo nivel de jusante da usina analisada
    if(tail(as.numeric(pats), 1) < limjus - 3) {

        # Caso limjus seja muito maior que os patamares, assume que todos sao sem remanso
        out <- rep(TRUE, length(pats))

        return(out)

    } else {

        # Do contrario, considera que todos os patamares menores ou no entorno de limjus sao sem
        # remanso de cara
        semrmns <- as.numeric(pats) < (limjus + 1)

        # Determina convergem no inicio do dominio
        semrmns <- v_convimed & semrmns
        semrmns[is.na(semrmns)] <- FALSE

        # Avalia a concavidade e determina quais nao tem influencia de remanso
        for(i in which(semrmns)) {

            # Pega dados e suavizacao do patamar qualificado em questao
            d <- dat$hist_est[pat == pats[i]]
            m <- tail(dat$patinfo[[i]]$tend, 1)[[1]]

            # Amostra valores no inicio do dominio e calcula segunda derivada
            v_x <- seq(d[valido == TRUE, min(vazao)], vazmin, length.out = 300)
            v_y <- predict(m, newdata = data.frame(vazao = v_x))
            v_d1 <- (v_y[-1]  - v_y[-length(v_y)]) / diff(v_x)[1]
            v_d2 <- (v_d1[-1] - v_d1[-length(v_d1)]) / diff(v_x)[1]

            # Determina se a concavidade e coerente com ausencia remanso
            if(mean(v_d2 < 0) > 0.9) semrmns[i] <- TRUE else semrmns[i] <- FALSE
        }

        # Assume que todos os patamares abaixo do maior sem influencia de remanso tambem nao sofrem influencia
        out <- cumsum(semrmns[seq(length(semrmns), 1)])
        out <- (out > 0)[seq(length(semrmns), 1)]

        return(out)
    }
}

tratavazconv <- function(dat) {

    patremns <- sapply(dat$patinfo, "[[", "remanso")

    out      <- sapply(dat$patinfo, "[[", "vazconv")
    vazconv <- out[patremns]

    # Corrige as convergencias para que sejam sempre crescentes
    for(i in seq(vazconv)) {

        if(is.na(vazconv[i])) next

        ultvazconv <- vazconv[1:(i - 1)]

        ultvazconv <- max(ultvazconv, na.rm = TRUE)

        # Testa se foi possivel achar tal vazao e, caso positivo, garante que aquela em questao e >=
        if((length(ultvazconv) != 0) && (vazconv[i] < ultvazconv)) vazconv[i] <- ultvazconv
    }

    out[patremns] <- vazconv

    return(out)
}
