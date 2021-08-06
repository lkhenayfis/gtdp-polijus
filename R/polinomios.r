################################ FUNCOES PARA AJUSTE DOS POLINOMIOS ################################

#' Ajuste de curvas polinomiais por partes
#' 
#' Genérica e métodos para ajuste das curvas polinomiais por partes a cada conjunto de dados
#' 
#' @param dat objeto do tipo \code{datcbase} ou \code{data.table}. Ver Detalhes
#' @param ext dado para extensão do ajuste. Ver Detalhes
#' @param graus vetor de até três posições indicando os graus de cada parte polinomial. Ver Detalhes
#' @param pto_turbmax vetor de duas posicoes indicando coordenadas da conexão entre polinômios 
#'     ajustados aos dados históricos
#' @param ... demais parâmetros específicos de cada método. Ver Detalhes
#' 
#' @return objeto \code{polijusU}
#' 
#' @seealso \code{polijusU} para detalhes do objeto retornado
#' 
#' @export

fitpoli <- function(dat, ext, graus, pto_turbmax, ...) UseMethod("fitpoli")

#' @param pto_ext vetor de duas posicoes indicando coordenadas da conexão entre último polinômio
#'     ajustados aos dados históricos e polinômio dos dados de extensão
#' 
#' @section Curva base:
#' 
#' Quando \code{dat} é da classe \code{datcbase} o ajuste da curva pode conter até três partes 
#' polinomiais de até quarto grau. Existem quatro possíveis configurações de ajuste
#' 
#' \itemize{
#'     \item Um polinômio para os dados históricos e nenhum de extensão (tipo \code{H1})
#'     \item Dois polinômios para os dados históricos e nenhum de extensão (tipo \code{H2})
#'     \item Um polinômio para os dados históricos e um para dados de extensão (tipo \code{H1E1})
#'     \item Dois polinômios para os dados históricos e um para dados de extensão (tipo \code{H2E1})
#' }
#' 
#' Qual ajuste será realizado é decidido internamente em função dos parâmetros passados para a 
#' função.
#' 
#' Ajustes do tipo \code{H1} exigem apenas \code{dat} e \code{graus}, sendo este um escalar (mais
#' posições serão ignoradas). Qualquer outro parâmetro fornecido além destes dois levará a função a 
#' tentar ajustar formas mais complexas
#' 
#' Ajustes do tipo \code{H2} demandam \code{dat}, \code{graus} e \code{pto_turbmax}. O dado de 
#' extensão deve ser fornecido ou como um inteiro ou string indicando qual dado da lista deve ser 
#' usado (use \code{polijus::plot.datcbase} para visualizar as opções). \code{graus} deve conter
#' dois elementos e \code{pto_turbmax} deve ser um vetor de duas posições indicando o ponto de 
#' conexão entre os dois polinômios (forma c(vazao, njus)).
#' 
#' Para ajustar uma configuração \code{H1E1}, os parâmetros necessários são \code{dat}, \code{ext} e
#' \code{graus} \code{pto_ext}. \code{ext} pode ser um inteiro ou strig indicando qual dos dados de 
#' extensão deve ser usado em conjunto com o histórico (use \code{polijus::plot.datcbase} para 
#' visualizar as opções). \code{graus} deve conter dois elementos e \code{pto_ext} indica o ponto de
#' conexão entre polinômio para histórico e polinômio de extensão (no mesmo molde 
#' \code{pto_turbmax}). 
#' 
#' Deve ser notado que o conjunto de parâmetros para \code{H2} e \code{H1E1} são similares bem 
#' similares. Para evitar problemas, caso seja fornecido um conjunto de parâmetros ambíguos à 
#' \code{fitpoli}, a função tentará realizar um ajuste do tipo \code{H1E1}, ignorando os parâmetros
#' associados ao tipo \code{H2}. Caso não haja parâmetros suficientes para realizar o \code{H1E1}, 
#' é retornado erro.
#' 
#' Finalmente, um ajuste do tipo \code{H2E1} demanda o fornecimento de todos os parâmetros 
#' anteriormente citados, da forma como foram detalhados, exceto por \code{graus} que agora deveria
#' conter 3 elementos.
#' 
#' @examples 
#' 
#' # METODO datcbase 
#' 
#' dbase <- extraibase(dummydata)
#' 
#' # ajuste de curva base H1
#' polibase <- fitpoli(dbase, graus = 3)
#' 
#' # ajuste de curva base H2
#' polibase <- fitpoli(dbase, graus = c(3, 3), pto_turbmax = c(400, 32.1))
#' 
#' # ajuste de curva base H1E1
#' polibase <- fitpoli(dat, ext = "CAD", graus = c(3, 3), pto_ext = c(1000, 34))
#' polibase <- fitpoli(dat, ext = 1, graus = c(3, 3), pto_ext = c(1000, 34)) # mesma coisa que acima
#' 
#' # ajuste de curva base H2E1
#' polibase <- fitpoli(dat, ext = "CAD", graus = c(3, 3, 3), pto_turbmax = c(400, 32.1), 
#'     pto_ext = c(1000, 34))
#' 
#' @rdname fitpoli
#' 
#' @export

fitpoli.datcbase <- function(dat, ext, graus, pto_turbmax, pto_ext) {

    l_parse <- parseargsbase(dat, ext, graus, pto_turbmax, pto_ext)
    l_func  <- l_parse[[1]]
    scales  <- l_parse[[2]]

    out <- eval(do.call(call, l_func))
    out <- rescale(out, scales)

    return(out)
}

# HELPERS ------------------------------------------------------------------------------------------

parseargsbase <- function(dat, ext, graus, pto_turbmax, pto_ext) {

    dat <- copybase(dat)

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

    dath1 <- dat$hist[vazao <= pto_turbmax[1]]

    vazrestrh1 <- seq(attr(dat, "vazzero"), pto_turbmax[1], length.out = 1000)

    Ah1 <- outer(dath1[, vazao], 0:graus[1], "^")
    bh1 <- data.matrix(dath1[, njus, drop = FALSE])
    Eh1 <- outer(tail(vazrestrh1, 1), 0:graus[1], function(x, y) x^y)
    fh1 <- pto_turbmax[2]
    Gh1 <- outer(vazrestrh1, 0:graus[1], function(x, y) y * x^(y - 1))
    hh1 <- matrix(rep(0, length(vazrestrh1)))

    dath2 <- dat$hist[vazao > pto_turbmax[1]]

    vazrestrh2 <- seq(pto_turbmax[1], max(dath2[, vazao]), length.out = 1000)

    Ah2 <- outer(dath2[, vazao], 0:graus[1], "^")
    bh2 <- data.matrix(dath2[, njus, drop = FALSE])
    Eh2 <- outer(head(vazrestrh2, 1), 0:graus[2], function(x, y) x^y)
    fh2 <- pto_turbmax[2]
    Gh2 <- outer(vazrestrh2, 0:graus[1], function(x, y) y * x^(y - 1))
    hh2 <- matrix(rep(0, length(vazrestrh2)))

    A <- Matrix::bdiag(Ah1, Ah2)
    b <- rbind(bh1, bh2)
    E <- Matrix::bdiag(Eh1, Eh2)
    f <- rbind(fh1, fh2)
    G <- Matrix::bdiag(Gh1, Gh2)
    h <- Matrix::Matrix(rbind(hh1, hh2))

    sol  <- limSolve::lsei(A, b, E, f, G, h, fulloutput = TRUE)

    breaks <- c(1, graus[1] + 1, graus[1] + 2, graus[1] + graus[2] + 2)
    coef   <- lapply(seq(graus), function(i) sol$X[breaks[2 * i - 1]:breaks[2 * i]])

    bounds <- list(range(vazrestrh1), range(vazrestrh2))

    new_polijusU(coef, bounds, dat, sol$covar, "H2", "curvabase")
}

fit_baseH1E1 <- function(dat, ext, graus, pto_ext, ...) {

    vazao <- njus <- NULL

    dath1 <- dat$hist

    vazrestrh1 <- seq(attr(dat, "vazzero"), pto_ext[1], length.out = 1000)

    Ah1 <- outer(dath1[, vazao], 0:graus[1], "^")
    bh1 <- data.matrix(dath1[, njus, drop = FALSE])
    Eh1 <- rbind(outer(tail(vazrestrh1, 1), 0:graus[1], function(x, y) x^y),
                outer(tail(vazrestrh1, 1), 0:graus[1], function(x, y) y * x^(y - 1)))
    fh1 <- rbind(pto_ext[2], 0)
    Gh1 <- outer(vazrestrh1, 0:graus[1], function(x, y) y * x^(y - 1))
    hh1 <- matrix(rep(0, length(vazrestrh1)))

    date1 <- dat$ext[[1]]

    vazrestrh2 <- seq(pto_ext[1], max(date1[, vazao]), length.out = 1000)

    Ae1 <- outer(date1[, vazao], 0:graus[1], "^")
    be1 <- data.matrix(date1[, njus, drop = FALSE])
    Ee1 <- rbind(outer(tail(vazrestrh2, 1), 0:graus[2], function(x, y) -y * x^(y - 1)),
                outer(head(vazrestrh2, 1), 0:graus[2], function(x, y) x^y))
    fe1 <- rbind(pto_ext[2])
    Ge1 <- outer(vazrestrh2, 0:graus[1], function(x, y) y * x^(y - 1))
    he1 <- matrix(rep(0, length(vazrestrh2)))

    A <- Matrix::bdiag(Ah1, Ae1)
    b <- rbind(bh1, be1)
    E <- cbind(rbind(Eh1, 0), rbind(0, Ee1))
    f <- rbind(fh1, fe1)
    G <- Matrix::bdiag(Gh1, Ge1)
    h <- Matrix::Matrix(rbind(hh1, he1))

    sol  <- limSolve::lsei(A, b, E, f, G, h, fulloutput = TRUE)

    breaks <- c(1, graus[1] + 1, graus[1] + 2, graus[1] + graus[2] + 2)
    coef   <- lapply(seq(graus), function(i) sol$X[breaks[2 * i - 1]:breaks[2 * i]])

    bounds <- list(range(vazrestrh1), range(vazrestrh2))

    new_polijusU(coef, bounds, dat, sol$covar, "H1E1", "curvabase")
}

fit_baseH2E1 <- function(dat, ext, graus, pto_turbmax, pto_ext) {

    vazao <- njus <- NULL

    dath1 <- dat$hist[vazao <= pto_turbmax[1]]

    vazrestrh1 <- seq(attr(dat, "vazzero"), pto_turbmax[1], length.out = 1000)

    Ah1 <- outer(dath1[, vazao], 0:graus[1], "^")
    bh1 <- data.matrix(dath1[, njus, drop = FALSE])
    Eh1 <- outer(tail(vazrestrh1, 1), 0:graus[1], function(x, y) x^y)
    fh1 <- pto_turbmax[2]
    Gh1 <- outer(vazrestrh1, 0:graus[1], function(x, y) y * x^(y - 1))
    hh1 <- matrix(rep(0, length(vazrestrh1)))

    dath2 <- dat$hist[vazao > pto_turbmax[1]]

    vazrestrh2 <- seq(pto_turbmax[1], pto_ext[1], length.out = 1000)

    Ah2 <- outer(dath2[, vazao], 0:graus[1], "^")
    bh2 <- data.matrix(dath2[, njus, drop = FALSE])
    Eh2 <- rbind(outer(head(vazrestrh2, 1), 0:graus[2], function(x, y) x^y),
                 outer(tail(vazrestrh2, 1), 0:graus[2], function(x, y) x^y),
                 outer(tail(vazrestrh2, 1), 0:graus[2], function(x, y) y * x^(y - 1)))
    fh2 <- rbind(pto_turbmax[2], pto_ext[2], 0)
    Gh2 <- outer(vazrestrh2, 0:graus[1], function(x, y) y * x^(y - 1))
    hh2 <- matrix(rep(0, length(vazrestrh2)))

    date1 <- dat$ext[[1]]

    vazrestre1 <- seq(pto_ext[1], max(date1[, vazao]), length.out = 1000)

    Ae1 <- outer(date1[, vazao], 0:graus[1], "^")
    be1 <- data.matrix(date1[, njus, drop = FALSE])
    Ee1 <- rbind(outer(tail(vazrestre1, 1), 0:graus[3], function(x, y) -y * x^(y - 1)),
                 outer(head(vazrestre1, 1), 0:graus[3], function(x, y) x^y))
    fe1 <- rbind(pto_ext[2])
    Ge1 <- outer(vazrestre1, 0:graus[1], function(x, y) y * x^(y - 1))
    he1 <- matrix(rep(0, length(vazrestre1)))

    A <- Matrix::bdiag(Ah1, Ah2, Ae1)
    b <- rbind(bh1, bh2, be1)
    E <- Matrix::bdiag(Eh1, cbind(rbind(Eh2, 0), rbind(0, 0, Ee1)))
    f <- rbind(fh1, fh2, fe1)
    G <- Matrix::bdiag(Gh1, Gh2, Ge1)
    h <- Matrix::Matrix(rbind(hh1, hh2, he1))

    sol  <- limSolve::lsei(A, b, E, f, G, h, fulloutput = TRUE)

    breaks <- c(1, graus[1] + 1,
                graus[1] + 2, graus[1] + graus[2] + 2,
                graus[1] + graus[2] + 3, graus[1] + graus[2] + graus[3] + 3)
    coef   <- lapply(seq(graus), function(i) sol$X[breaks[2 * i - 1]:breaks[2 * i]])

    bounds <- list(range(vazrestrh1), range(vazrestrh2), range(vazrestre1))

    new_polijusU(coef, bounds, dat, sol$covar, "H1E1", "curvabase")
}