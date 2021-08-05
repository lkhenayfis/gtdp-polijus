
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

#' @export 

print.polijusU <- function(polijus, ...) {

    out <- coef(polijus)
    cli::cli_h2(attr(polijus, "tag"))
    print(out)
}

#' @export 

coef.polijusU <- function(polijus, ...) {

    npoli <- attr(polijus, "npoli")
    out <- lapply(seq(npoli), function(i) {
        dom <- polijus$bounds[[i]]
        coef <- polijus$coefs[[i]]

        matrix(c(dom, coef), nrow = 1,
            dimnames = list(paste0("poli", i), c("Vaz_min", "Vaz_max", paste0("A", 0:4))))
    })

    out <- do.call(rbind, out)

    return(out)
}

#' @export 

fitted.polijusU <- function(polijus, ...) {

    vazao <- poli <- NULL

    npoli <- attr(polijus, "npoli")

    coefs <- polijus$coefs

    fitted <- lapply(seq(npoli), function(i) {
        vaz    <- polijus$model[poli == i, vazao]
        coefI  <- coefs[[i]]
        fit    <- lapply(seq(coefI), function(i) vaz^(i - 1) * coefI[i])

        rowSums(do.call(cbind, fit))
    })

    fitted <- unlist(fitted)

    return(fitted)
}

#' @export 

predict.polijusU <- function(polijus, newdata, ...) {

    vazao <- NULL

    setDT(newdata)

    npoli <- attr(polijus, "npoli")

    coefs <- polijus$coefs

    predicted <- lapply(seq(npoli), function(i) {
        bounds <- polijus$bounds[[i]]
        vaz    <- newdata[vazao %between% bounds, vazao]
        coefI  <- coefs[[i]]
        pred   <- lapply(seq(coefI), function(i) vaz^(i - 1) * coefI[i])

        rowSums(do.call(cbind, pred))
    })

    predicted <- unlist(predicted)
    predicted <- predicted[!duplicated(predicted)] # caso haja pontos exatamente iguais aos bounds (pouco provavel)

    return(predicted)
}

#' @export

rescale <- function(polijus, scales, inv = FALSE) {

    npoli <- attr(polijus, "npoli")

    coefs <- polijus$coefs

    mu  <- scales[[1]]
    sig <- scales[[2]]

    if(!inv) {

        DELTA <- lapply(seq(npoli), function(i) {
            coefI <- coefs[[i]]
            outer(1:5, 1:5, function(m, n) {
                (n >= m) * sig[2] * (-1)^(n - m) * choose(n - 1, m - 1) * mu[1]^(n - m) / sig[1]^(n - 1) * coefI[n]
            })
        })
        d <- lapply(seq(npoli), function(i) c(mu[2], double(4)))

        polijus$coefs  <- lapply(seq(npoli), function(i) c(DELTA[[i]] %*% coefs[[i]] + d[[i]]))
        polijus$bounds <- lapply(seq(npoli), function(i) polijus$bounds[[i]] * sig[1] + mu[1])
        polijus$model[, 1:2 := mapply(.SD, mu, sig, FUN = function(d, u, s) d * s + u,
            SIMPLIFY = FALSE), .SDcols = 1:2]

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
                (n >= m) * choose(n - 1, m - 1) * mu[1]^(n - m) * sig[1]^(n - 1) * coefI[n]
            })
        })
        d <- lapply(seq(npoli), function(i) c(-mu[2] / sig[2], double(4)))

        polijus$coefs  <- lapply(seq(npoli), function(i) c(DELTA[[i]] %*% coefs[[i]] + d[[i]]) / sig[2])
        polijus$bounds <- lapply(seq(npoli), function(i) (polijus$bounds[[i]] - mu[1]) / sig[1])
        polijus$model[, 1:2 := mapply(.SD, mu, sig, FUN = function(d, u, s) (d - u) / s,
            SIMPLIFY = FALSE), .SDcols = 1:2]

        bDELTA  <- as.matrix(do.call(Matrix::bdiag, DELTA))
        bDELTAt <- as.matrix(do.call(Matrix::bdiag, lapply(DELTA, t)))
        dn <- dimnames(polijus$vcov)
        polijus$vcov <- bDELTA %*% polijus$vcov %*% bDELTAt
        dimnames(polijus$vcov) <- dn

        polijus$fitted <- fitted(polijus)
    }

    return(polijus)
}
