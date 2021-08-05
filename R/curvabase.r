########################## FUNCOES PARA EXTRACAO E EXTENSAO DA CURVA BASE ##########################

#' Extrator de curva abase
#' 
#' Facilita a extração dos pontos relativos à curva base a partir do histórico estável
#' 
#' Por padrão a curva base é composta pelos dados sem influência de remanso como identificado
#' anteriormente em \code{\link{evalremanso}}. Em determinados casos de dados menos consistentes 
#' este processo de identificação não resulta em um conjunto representativo para a curva base, sendo
#' necessário intervir manualmente.
#' 
#' Para este propósito existe o parâmetro \code{subset}. Se for fornecido, a identificação padrão de
#' dados da curva base é completamente ignorada e são utilizados apenas os registros informados
#' através de \code{subset}.
#' 
#' @param dat objeto \code{datpoli} do qual extrair dados para base
#' @param subset opcional -- vetor de inteiros ou lógico indicando as posições do dado estável para
#'     usar na curva base. Ver Detalhes
#' 
#' @return objeto \code{datcbase} -- lista de dois elementos:
#'     \itemize{
#'         \item no primeiro está a parte do dado histórico utilizada para curva base
#'         \item no segundo uma lista de alternativas de dados para extensão do ajuste
#' }
#' 
#' @export

extraibase <- function(dat, subset = NULL) {

    vazao <- valido <- pat <- base <- NULL

    if(!missing("subset")) {
        dbase <- dat$hist_est[subset]
    } else {
        remanso <- sapply(dat$patinfo, "[[", "remanso")
        vazconv <- sapply(dat$patinfo, "[[", "vazconv")

        dbase <- copy(dat$hist_est)

        # Dados sem efeito de remanso
        dbase[, base := FALSE]
        dbase[pat %in% names(remanso[!remanso]), base := TRUE]
        for(p in names(remanso[remanso])) {
            dbase[(pat == p) & (vazao > vazconv[p]), base := TRUE]
        }

        dbase <- dbase[(valido == TRUE) & (base == TRUE)]
    }

    maxvaz    <- max(3 * dbase$vazao, attr(dat, "vazmax"))
    vazextrap <- seq(max(dbase$vazao), maxvaz, length.out = 500)
    lextrap   <- lapply(dat$ext, function(f) {
        data.table(vazao = vazextrap, njus = f(vazextrap))
    })

    dbase <- list(hist = dbase, ext = lextrap)
    class(dbase) <- "datcbase"
    attr(dbase, "vazzero") <- 0

    return(dbase)
}

#' Extrapolação logarítmica
#' 
#' @export

extraplog <- function(dbase, tol = 1e-5) {

    if(!class(dbase) == "datcbase") {
        stop(paste0("dbase fornecido nao tem classe 'datcbase' -- ",
            "Use polijus::extrabase para obter o dado corretamente"))
    }

    dextrap <- copy(dbase[[1]])
    setorder(dextrap, vazao)
    dextrap <- dextrap[vazao > median(vazao)]
    dextrap[, c("vazao", "njus") := list(log(vazao), log(njus))]

    vazhistmin <- min(log(dbase[[1]]$vazao))

    vazextrap <- dbase$ext[[1]]$vazao

    vazteste <- seq(dextrap$vazao[1], max(dextrap$vazao), length.out = 500)
    for(vazi in vazteste) {

        if(dextrap[vazao > vazi, .N == 0]) stop("Nao foi possivel extrapolar o dado. Tente aumentar a tolerancia")

        pesos <- dextrap[vazao > vazi, (vazao - vazhistmin) / (max(vazao) - vazhistmin)]

        reg <- lm(njus ~ vazao, dextrap[vazao > vazi], weights = pesos)

        erro <- sum(reg$residuals^2)

        if(erro < tol) break
    }

    plot(log(dbase$hist$vazao), log(dbase$hist$njus), col = "deepskyblue2", pch = 16,
        panel.first = grid(col = "grey85"),
        xlab = expression("Log Vazão [m"^3 * "/s]"), ylab = "Log Nível de jusante [m]",
        main = "Região linear extrapolada")
    points(dextrap[vazao > vazi, c("vazao", "njus")], col = "purple3", pch = 16)
    abline(reg, lwd = 2)
    legend("bottomright", inset = 0.02,
            legend = c("Dados históricos", "Região linear", "Ajuste linear"),
            pch = c(16, 16, NA),
            lty = c(NA, NA, 1),
            lwd = c(NA, NA, 2),
            col = c("deepskyblue2", "purple3", "black"))

    vazao <- log(vazextrap)
    njus  <- predict(reg, newdata = data.frame(vazao = vazao))
    out   <- data.table(vazao = exp(vazao), njus = exp(njus))

    dbase[[2]]$EXTRAP <- out

    return(dbase)
}

# HELPERS ------------------------------------------------------------------------------------------

#' @export

copy.datcurvavase <- function(dat) {
    out <- list(hist = NULL, ext = NULL)
    out$hist <- copy(dat$hist)
    out$ext  <- lapply(dat$ext, copy)
    return(out)
}

#' @export

scale.datcbase <- function(dat, center = TRUE, scale = TRUE) {

    if(center) {
        medias <- sapply(c("vazao", "njus"), function(col) {
            vhist <- dat[[1]][[col]]
            lext  <- lapply(dat[[2]], function(ext) ext[[col]])
            mean(c(vhist, unlist(lext)))
        })

        dat$hist[, 3:2 := mapply(.SD, medias, SIMPLIFY = FALSE, FUN = function(v, m) v - m), .SDcols = 3:2]
        lapply(dat$ext, function(d) {
            d[, 1:2 := mapply(.SD, medias, SIMPLIFY = FALSE, FUN = function(v, m) v - m), .SDcols = 1:2]
            invisible(NULL)
        })
    }

    if(scale) {
        sds <- sapply(c("vazao", "njus"), function(col) {
            vhist <- dat[[1]][[col]]
            lext  <- lapply(dat[[2]], function(ext) ext[[col]])
            sd(c(vhist, unlist(lext)))
        })

        dat$hist[, 3:2 := mapply(.SD, sds, SIMPLIFY = FALSE, FUN = function(v, m) v / m), .SDcols = 3:2]
        lapply(dat$ext, function(d) {
            d[, 1:2 := mapply(.SD, sds, SIMPLIFY = FALSE, FUN = function(v, m) v / m), .SDcols = 1:2]
            invisible(NULL)
        })
    }

    return(list(medias, sds))
}
