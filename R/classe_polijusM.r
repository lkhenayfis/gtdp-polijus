############################### CLASSE POLIJUS MULTIVARIADO E METODOS ##############################

#' Objeto \code{polijusM}
#' 
#' Contrutor e métodos da classe \code{polijusM} -- curva de jusante multivariada
#' 
#' @name polijusM
NULL

#' @details Construtor primário de objetos \code{polijusM}
#' 
#' @param datorig objeto \code{datpoli} contendo os dados originais utilizados para ajuste
#' @param polibase curva base ajustada
#' @param l_poliind lista de curvas individuais ajustadas
#' 
#' @examples 
#' 
#' # ajusta uma curva base --------------------------------
#' datbase <- extraibase(dummydata)
#' 
#' ext   <- 1
#' graus <- list(2:4, 2:4)
#' 
#' maxvaz   <- 1.1 * dbase$hist[, max(vazao)]
#' pto_ext0 <- c(maxvaz, dummydata$ext[[ext]](maxvaz))
#' 
#' polibase <- optpoli(datbase, ext, graus, pto_ext0 = pto_ext0)[[2]]
#' 
#' # seleciona e ajusta curvas individuais ----------------
#' 
#' l_datind <- extraipats(dummydata, polibase, 10)
#' l_poliind <- lapply(l_datind, function(dat) {
#'     optpoli(dat, polibase, list(2:4), vaz_ext0 = dat[, max(vazao)])
#' })
#' 
#' # monta objeto final -----------------------------------
#' 
#' polijus <- new_polijusM(dummydata, polibase, l_poliind)
#' 
#' \dontrun{
#'     print(polijus)
#'     coef(polijus)
#'     summary(polijus)
#'     plot(polijus)
#' }
#' 
#' @return objeto \code{polijusM}, uma lista contendo
#' \describe{
#'     \item{curvas}{lista de curvas ajustadas ao dado}
#'     \item{model}{dado de vazao, nível de jusante e monstante representados pela função ajustada}
#' }
#' 
#' Adicionalmente possui atributos
#' \describe{
#'     \item{cod}{código da usina ajustada}
#'     \item{ncurvas}{número de curvas ajustadas}
#'     \item{npolis}{número de polinômios por curva}
#'     \item{refs}{níveis de montante de referência para cada curva ajustda}
#' }
#' 
#' @rdname polijusM
#' 
#' @export

new_polijusM <- function(datorig, polibase, l_poliind) {

    vazao <- njus <- nmont <- NULL

    hist <- datorig$hist_est[valido == TRUE, .(vazao, njus, nmont)]

    curvas  <- c(list(curvabase = polibase), l_poliind)
    ncurvas <- length(curvas)
    npolis  <- sapply(curvas, function(poli) nrow(coef(poli)))

    refind <- sapply(l_poliind, function(poli) as.numeric(sub(".* ", "", attr(poli, "tag"))))

    desloc  <- mapply(l_poliind, refind, FUN = function(poli, ref) ref - coef(poli)[1, 3])
    refbase <- coef(polibase)[1, 3] + mean(desloc)

    refs <- structure(c(refbase, refind), names = names(npolis))

    out <- list(curvas = curvas, model = hist)

    attr(out, "cod")     <- attr(datorig, "cod")
    attr(out, "ncurvas") <- ncurvas
    attr(out, "npolis")  <- npolis
    attr(out, "refs")    <- refs
    class(out) <- "polijusM"

    return(out)
}

# METODOS ------------------------------------------------------------------------------------------

#' @param x objeto \code{polijusM}
#' @param ... demais parametros
#' 
#' @rdname polijusM
#' 
#' @export

print.polijusM <- function(x, ...) {

    cod <- attr(x, "cod")
    cli::cli_h2(c("Usina ", cod))

    curvas <- unname(sapply(x$curvas, attr, "tag"))
    npoli  <- unname(attr(x, "npolis"))
    refs   <- unname(attr(x, "refs"))
    mat    <- data.frame(Curva = curvas, Polinomios = npoli, Referencia = refs, row.names = NULL)

    print(mat)
}

#' @rdname polijusM
#' 
#' @export

coef.polijusM <- function(x, ...) {

    coefs <- lapply(x$curvas, coef)
    for(i in seq(coefs)) {
        mat <- coefs[[i]]
        rownames(mat) <- NULL

        curva <- rep(i, nrow(mat))
        poli  <- seq(nrow(mat))

        pre <- matrix(c(curva, poli), nrow(mat), dimnames = list(NULL, c("Curva", "Polinomio")))
        coefs[[i]] <- cbind(pre, mat)
    }
    coefs <- do.call(rbind, coefs)

    print(coefs)
}

#' @param object objeto tipo \code{polijusM}
#' @param ... demais parâmetros 
#' 
#' @rdname polijusM
#' 
#' @export 

fitted.polijusM <- function(object, ...) {

    model   <- copy(object$model)
    ncurvas <- attr(object, "ncurvas")
    refs    <- attr(object, "refs")

    intervalo <- findInterval(model$nmont, refs)

    curvainf <- intervalo
    curvainf[curvainf == 0] <- 1

    curvasup <- intervalo + 1
    curvasup[curvasup > ncurvas] <- ncurvas

    predcurva <- function(x, ind) predict(object$curvas[[ind[1]]], newdata = data.frame(vazao = x))

    model[, c("cinf", "csup") := .(curvainf, curvasup)]
    model[, c("refinf", "refsup") := .(refs[cinf], refs[csup])]
    model[, predinf := predcurva(vazao, cinf), by = cinf]
    model[, predsup := predcurva(vazao, csup), by = csup]

    predfinal <- function(x, x1, x2, y1, y2) {
        b1 <- (y2 - y1) / (x2 - x1)
        ifelse(is.na(b1), y1, y1 + b1 * (x - x1))

    }

    model[, pred := predfinal(nmont, refinf, refsup, predinf, predsup)]

    return(model$pred)
}

#' @param newdata data.frame ou data.table opcional contendo dados com os quais realizar previsão
#' 
#' @rdname polijusM
#' 
#' @export 

predicted.polijusM <- function(object, newdata, ...) {

    object$model <- newdata
    fitted(object)
}

#' @rdname polijusM
#' 
#' @export

residuals.polijusM <- function(object, ...) {

    fitt <- fitted(object)
    object$model$njus - fitt
}

#' @rdname polijusM
#' 
#' @export

summary.polijusM <- function(object, ...) {

    curvas <- object$curvas
    names(curvas)[1] <- "object"
    do.call(summary, curvas)
}