################################ FUNCOES PARA AJUSTE DOS POLINOMIOS ################################

#' @export

fitpoli <- function(dat, ext, graus, pto_turbmax, pto_ext, opcoes) UseMethod("fitpoli")

#' @export

fitpoli.datcbase <- function(dat, ext, graus, pto_turbmax, pto_ext, opcoes) {

    l_parse <- parseargsbase(dat, ext, graus, pto_turbmax, pto_ext)
    call    <- l_parse[[1]]
    scales  <- l_parse[[2]]

}

# HELPERS ------------------------------------------------------------------------------------------

parseargsbase <- function(dat, ext, graus, pto_turbmax, pto_ext) {

    dat <- copy(dat)

    ngraus <- length(graus)
    if(ngraus > 3) {
        warning("'graus' possui mais de tres elementos -- reduzindo para apenas os tres primeiros")
        graus  <- graus[1:3]
        ngraus <- 3
    }

    temext <- !missing("ext")
    temptm <- !missing("pto_turbmax")
    tempve <- !missing("pto_ext")

    if(!(tempve == temext)) {
        warning("Foi fornecido apenas um de 'pto_ext' e 'ext' --- ",
            "\nSe deseja ajustar um polinomio para dados de extensao, selecione o conjunto atraves",
            " do parametro 'ext' e um ponto de conexao com os dados historicos atraves",
            " do parametro 'pto_ext' (veja ?polijus::fitpoli)")
        temext <- FALSE
        tempve <- FALSE
    }

    # Identificacao de tipo de ajuste
    if((ngraus == 1)) {

        # Caso de polinomio unico - ignora todos os parametros extras
        extras <- structure(c(temext, temptm, tempve), names = c("ext", "pto_turbmax", "pto_ext"))
        if(any(extras)) {
            ignore <- names(extras)[extras]
            warning("'graus' contem apenas um elemento, indicando polinomio unico -- os parametros (",
                ignore, ") serao ignorados", "\nPara ajuste de multiplos polinomios forneca ate tres ",
                "elementos em graus, acompanhados dos parametros adequados (veja ?polijus::fitpoli)")
        }

        func <- "fit_baseH1"
        ext  <- NULL
        pto_turbmax <- c(NA_real_, NA_real_)
        pto_ext  <- c(NA_real_, NA_real_)

    } else if(ngraus == 2) {

        # Caso de dois polinomios -- pode ser H1E1 ou H2, dando preferencia para o primeiro

        if(temext) {

            if(temptm) {
                warning("'graus' contem 2 elementos, indicando ajuste de dois polinomios, porem tanto ",
                    "'ext' quanto 'pto_turbmax' foram fornecidos -- 'pto_turbmax' sera ignorado em favor",
                    " do ajuste de um polinomio para dados historicos e outro para dados de extensao",
                    "\n Se deseja ajustar dois polinomios para dados historicos e nenhum para exntensao,",
                    " forneca apenas 'pto_turbmax' (veja ?polijus::fitpoli)")
            }

            func <- "fit_baseH1E1"
            pto_turbmax <- c(NA_real_, NA_real_)

        } else {

            if(!temptm) {
                stop("'graus' contem 2 elementos, indicando ajuste de dois polinomios, porem nem ",
                    "'ext' nem 'pto_turbmax' foram fornecidos",
                    "\n Se deseja ajustar dois polinomios para dados historicos ou um para historico ",
                    "e outro para extensao, forneca 'pto_turbmax' ou 'ext', respectivamente (veja ?polijus::fitpoli)")
            }

            func <- "fit_baseH2"
            ext  <- NULL
            pto_ext <- c(NA_real_, NA_real_)
        }
    } else {

        # Caso de tres polinomios

        if(!temext | !temptm) {
            stop("'graus' contem 3 elementos, indicando ajuste de tres polinomios, porem ",
                "'ext' e 'pto_turbmax' nao foram fornecidos",
                "\n Se deseja ajustar dois polinomios para dados historicos ",
                "e outro para extensao, forneca 'pto_turbmax' e 'ext' (veja ?polijus::fitpoli)")
        }

        func <- "fit_baseH2E1"
    }

    dat$ext <- dat$ext[ext]

    # OPERA POR REFERENCIA EM DAT E RETORNA LISTA DE PARAMTROS DE ESCALONAMENTO
    scales <- scale(dat)

    pto_turbmax <- (pto_turbmax - scales[[1]]) / scales[[2]]
    pto_ext     <- (pto_ext - scales[[1]]) / scales[[2]]

    attr(dat, "vazzero") <- -scales[[1]][1] / scales[[2]][1]

    call <- list(func, dat = dat, ext = ext, graus = graus, pto_turbmax = pto_turbmax, pto_ext = pto_ext)

    return(list(call, scales))
}

fit_baseH1 <- function(dat, graus, ...) {

    vazrestr <- seq(attr(dat, "vazzero"), max(dat$hist[, vazao]), length.out = 1000)

    A <- outer(dat$hist[, vazao], 0:graus[1], "^")
    b <- data.matrix(dat$hist[, njus, drop = FALSE])
    E <- NULL
    f <- NULL
    G <- outer(vazrestr, 0:graus[1], function(x, y) y * x^(y - 1))
    h <- matrix(rep(0, length(vazrestr)))

    coef   <- limSolve::lsei(A, b, E, f, G, h, fulloutput = TRUE)
    bounds <- dat$hist[, range(vazao)]

    new_polijusU(coef$X, bounds, dat, coef$covar, "curva base")
}