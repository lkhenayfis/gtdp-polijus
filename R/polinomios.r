################################ FUNCOES PARA AJUSTE DOS POLINOMIOS ################################

#' @export

fitpoli <- function(dat, ext, graus, pto_turbmax, pto_ext, opcoes) UseMethod("fitpoli")

#' @export

fitpoli.datcbase <- function(dat, ext, graus, pto_turbmax, pto_ext, opcoes) {

    l_parse <- parseargsbase(dat, ext, graus, pto_turbmax, pto_ext)
    l_func  <- l_parse[[1]]
#    scales  <- l_parse[[2]]

    out <- eval(do.call(call, l_func))
#    out <- rescale(out, scales)

    return(out)
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

    vazao <- njus <- NULL

    vazrestr <- seq(attr(dat, "vazzero"), max(dat$hist[, vazao]), length.out = 1000)

    A <- outer(dat$hist[, vazao], 0:graus[1], "^")
    b <- data.matrix(dat$hist[, njus, drop = FALSE])
    E <- NULL
    f <- NULL
    G <- outer(vazrestr, 0:graus[1], function(x, y) y * x^(y - 1))
    h <- matrix(rep(0, length(vazrestr)))

    sol    <- limSolve::lsei(A, b, E, f, G, h, fulloutput = TRUE)
    bounds <- range(vazrestr)

    new_polijusU(sol$X, bounds, dat, sol$covar, "H1", "curva base")
}

fit_baseH2 <- function(dat, graus, pto_turbmax, ...) {

    vazao <- njus <- NULL

    dat1 <- dat$hist[vazao <= pto_turbmax[1]]

    vazrestr1 <- seq(attr(dat, "vazzero"), max(dat1[, vazao]), length.out = 1000)

    A1 <- outer(dat1[, vazao], 0:graus[1], "^")
    b1 <- data.matrix(dat1[, njus, drop = FALSE])
    E1 <- outer(tail(vazrestr1, 1), 0:graus[1], function(x, y) x^y)
    f1 <- pto_turbmax[2]
    G1 <- outer(vazrestr1, 0:graus[1], function(x, y) y * x^(y - 1))
    h1 <- matrix(rep(0, length(vazrestr1)))

    dat2 <- dat$hist[vazao > pto_turbmax[1]]

    vazrestr2 <- seq(max(dat1[, vazao]), max(dat2[, vazao]), length.out = 1000)

    A2 <- outer(dat2[, vazao], 0:graus[1], "^")
    b2 <- data.matrix(dat2[, njus, drop = FALSE])
    E2 <- outer(head(vazrestr2, 1), 0:graus[2], function(x, y) x^y)
    f2 <- pto_turbmax[2]
    G2 <- outer(vazrestr2, 0:graus[1], function(x, y) y * x^(y - 1))
    h2 <- matrix(rep(0, length(vazrestr2)))

    A <- Matrix::bdiag(A1, A2)
    b <- rbind(b1, b2)
    E <- Matrix::bdiag(E1, E2)
    f <- rbind(f1, f2)
    G <- Matrix::bdiag(G1, G2)
    h <- Matrix::Matrix(rbind(h1, h2))

    sol  <- limSolve::lsei(A, b, E, f, G, h, fulloutput = TRUE)

    breaks <- c(1, graus[1] + 1, graus[1] + 2, graus[1] + graus[2] + 2)
    coef   <- lapply(seq(graus), function(i) sol$X[breaks[2 * i - 1]:breaks[2 * i]])

    bounds <- list(range(vazrestr1), range(vazrestr2))

    new_polijusU(coef, bounds, dat, sol$covar, "H2", "curvabase")
}
