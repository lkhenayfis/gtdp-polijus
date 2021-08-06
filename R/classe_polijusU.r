################################ CLASSE POLIJUS UNIVARIADO E METODOS ###############################

#' Objeto \code{polijusU}
#' 
#' Contrutor e métodos da classe \code{polijusU} -- curva de jusante univariada
#' 
#' @name polijusU
NULL

#' @details \code{new_polijusU} não deveria ser utilizada diretamente pelo usuário
#' 
#' @param coefs lista de um, dois ou três elementos contendo coeficientes ajustados das partes
#'     polinomiais
#' @param bounds lista contendo vetores de duas posições indicando o domínio de cada parte
#' @param dat objeto passado à \code{fitpoli}
#' @param vcov matriz de variâncias e covariâncas dos coeficientes ajustados
#' @param tipo o tipo de ajuste
#' @param tag uma tag para identificação de qual curva se trata no print
#' 
#' @return objeto \code{polijusU}, uma lista contendo
#' \describe{
#'     \item{coef}{lista de coeficientes por parte polinomial ajustada}
#'     \item{bounds}{lista de limites do domínio de cada parte polinomial ajustada}
#'     \item{model}{dado de vazao e nível ajustado com coluna de indicador da parte polinomial}
#'     \item{vcov}{matriz de variâncias e covariâncias dos coeficientes ajustados}
#'     \item{fitted}{vetor de pontos ajustados}
#' }
#' 
#' Adicionalmente possui atributos
#' \describe{
#'     \item{npoli}{número de partes polinomiais ajustadas}
#'     \item{graus}{vetor de graus das partes polinomiais}
#'     \item{tipo}{o tipo de ajuste realizado}
#'     \item{tag}{string indicando o nome do ajuste (curva base, patamar XXX, etc)}
#' }
#' 
#' @rdname polijusU
#' 
#' @export

new_polijusU <- function(coefs, bounds, dat, vcov, tipo, tag) {

    vazao <- poli <- NULL

    if(!is.list(coefs)) {
        coefs <- list(poli1 = coefs)
    } else {
        names(coefs) <- paste0("poli", seq(length(coefs)))
    }
    graus <- sapply(coefs, function(coef) length(coef) - 1)
    coefs <- lapply(coefs, function(coef) c(coef, double(5 - length(coef))))

    npoli <- length(coefs)

    if(!is.list(bounds)) {
        bounds <- list(poli1 = bounds)
    } else {
        names(bounds) <- paste0("poli", seq(length(bounds)))
    }

    dat$hist <- dat$hist[, c("vazao", "njus")]
    if(length(dat$ext) > 0) dat$ext <- dat$ext[[1]] else dat$ext <- NULL
    model  <- rbindlist(dat)
    breaks <- unlist(bounds)
    model[, poli := findInterval(vazao, breaks[!duplicated(breaks)], all.inside = TRUE)]
    setorder(model, vazao)

    if(is.null(vcov)) vcov <- matrix(0, graus[1] + 1, graus[1] + 1)
    nomes <- unlist(lapply(seq(npoli), function(i) paste0(i, "_", "A", 0:graus[i])))
    dimnames(vcov) <- list(nomes, nomes)

    nomes2 <- unlist(lapply(seq(npoli), function(i) paste0(i, "_", "A", 0:4)))
    aux <- matrix(.0, 5 * npoli, 5 * npoli, dimnames = list(nomes2, nomes2))
    aux[nomes, nomes] <- vcov

    out <- list(coefs = coefs, bounds = bounds, model = model, vcov = aux)
    attr(out, "npoli") <- npoli

    fitted <- fitted.polijusU(out)

    out$fitted <- fitted

    attr(out, "graus") <- graus
    attr(out, "tipo")  <- tipo
    attr(out, "tag")   <- tag
    class(out) <- "polijusU"

    return(out)
}

# METODOS ------------------------------------------------------------------------------------------

#' @param object,x objeto tipo \code{polijusU}
#' @param ... demais parâmetros 
#' 
#' @rdname polijusU
#' 
#' @export 

print.polijusU <- function(x, ...) {

    out <- coef(x)
    cli::cli_h2(attr(x, "tag"))
    print(out)
}

#' @rdname polijusU
#' 
#' @export 

coef.polijusU <- function(object, ...) {

    npoli <- attr(object, "npoli")
    out <- lapply(seq(npoli), function(i) {
        dom <- object$bounds[[i]]
        coef <- object$coefs[[i]]

        matrix(c(dom, coef), nrow = 1,
            dimnames = list(paste0("poli", i), c("Vaz_min", "Vaz_max", paste0("A", 0:4))))
    })

    out <- do.call(rbind, out)

    return(out)
}

#' @rdname polijusU
#' 
#' @export 

fitted.polijusU <- function(object, ...) {

    vazao <- poli <- NULL

    npoli <- attr(object, "npoli")

    coefs <- object$coefs

    fitted <- lapply(seq(npoli), function(i) {
        vaz    <- object$model[poli == i, vazao]
        coefI  <- coefs[[i]]
        fit    <- lapply(seq(coefI), function(i) vaz^(i - 1) * coefI[i])

        rowSums(do.call(cbind, fit))
    })

    fitted <- unlist(fitted)

    return(fitted)
}

#' @param newdata data.frame ou data.table opcional contendo dados com os quais realizar previsão
#' 
#' @rdname polijusU
#' 
#' @export 

predict.polijusU <- function(object, newdata, ...) {

    vazao <- NULL

    npoli <- attr(object, "npoli")
    coefs <- object$coefs
    bounds <- object$bounds

    setDT(newdata)
    breaks <- unlist(bounds)
    newdata[, poli := findInterval(vazao, breaks[!duplicated(breaks)], all.inside = TRUE)]

    predicted <- lapply(seq(npoli), function(i) {
        vaz    <- newdata[poli == i, vazao]
        coefI  <- coefs[[i]]
        pred   <- lapply(seq(coefI), function(i) vaz^(i - 1) * coefI[i])

        rowSums(do.call(cbind, pred))
    })

    predicted <- unlist(predicted)
    predicted <- predicted[!duplicated(predicted)] # caso haja pontos exatamente iguais aos bounds (pouco provavel)

    return(predicted)
}

# HELPERS ------------------------------------------------------------------------------------------

rescale <- function(polijus, scales, inv = FALSE) {

    npoli <- attr(polijus, "npoli")

    coefs <- polijus$coefs

    mu  <- scales[[1]]
    sig <- scales[[2]]

    if(!inv) {

        DELTA <- lapply(seq(npoli), function(i) {
            coefI <- coefs[[i]]
            outer(1:5, 1:5, function(m, n) {
                (n >= m) * sig[2] * (-1)^(n - m) * choose(n - 1, m - 1) * mu[1]^(n - m) / sig[1]^(n - 1)
            })
        })
        d <- lapply(seq(npoli), function(i) c(mu[2], double(4)))

        polijus$coefs  <- lapply(seq(npoli), function(i) c(DELTA[[i]] %*% coefs[[i]] + d[[i]]))
        polijus$bounds <- lapply(seq(npoli), function(i) polijus$bounds[[i]] * sig[1] + mu[1])

        # nao usa modificacao in-place de proposito para evitar modificar o argumento passado por
        # referencia no parent frame
        polijus$model[, c("vazao", "njus")] <- mapply(polijus$model[, c("vazao", "njus")], mu, sig,
            SIMPLIFY = FALSE, FUN = function(d, u, s) d * s + u)

        bDELTA  <- as.matrix(do.call(Matrix::bdiag, DELTA))
        bDELTAt <- as.matrix(do.call(Matrix::bdiag, lapply(DELTA, t)))
        dn <- dimnames(polijus$vcov)
        polijus$vcov <- bDELTA %*% polijus$vcov %*% bDELTAt
        dimnames(polijus$vcov) <- dn

        polijus$fitted <- fitted(polijus)

    } else {

        DELTA <- lapply(seq(npoli), function(i) {
            coefI <- coefs[[i]]
            outer(1:5, 1:5, function(m, n) {
                (n >= m) * sig[2]^(-1) * choose(n - 1, m - 1) * mu[1]^(n - m) * sig[1]^(m - 1)
            })
        })
        d <- lapply(seq(npoli), function(i) c(-mu[2] / sig[2], double(4)))

        polijus$coefs  <- lapply(seq(npoli), function(i) c(DELTA[[i]] %*% coefs[[i]] + d[[i]]))
        polijus$bounds <- lapply(seq(npoli), function(i) (polijus$bounds[[i]] - mu[1]) / sig[1])

        # nao usa modificacao in-place de proposito para evitar modificar o argumento passado por
        # referencia no parent frame
        polijus$model[, c("vazao", "njus")] <- mapply(polijus$model[, c("vazao", "njus")], mu, sig,
            SIMPLIFY = FALSE, FUN = function(d, u, s) (d - u) / s)

        bDELTA  <- as.matrix(do.call(Matrix::bdiag, DELTA))
        bDELTAt <- as.matrix(do.call(Matrix::bdiag, lapply(DELTA, t)))
        dn <- dimnames(polijus$vcov)
        polijus$vcov <- bDELTA %*% polijus$vcov %*% bDELTAt
        dimnames(polijus$vcov) <- dn

        polijus$fitted <- fitted(polijus)
    }

    return(polijus)
}
