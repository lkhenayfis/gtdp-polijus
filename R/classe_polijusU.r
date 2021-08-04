
#' @export

new_polijusU <- function(coefs, bounds, dat, vcov, tag) {

    if(!is.list(coefs)) {
        coefs <- list(poli1 = coefs)
    } else {
        names(coefs) <- paste0("poli", seq(length(coefs)))
    }
    coefs <- lapply(coefs, function(coef) c(coef, double(5 - length(coef))))

    if(!is.list(bounds)) {
        bounds <- list(poli1 = bounds)
    } else {
        names(bounds) <- paste0("poli", seq(length(bounds)))
    }

    dat$hist <- dat$hist[, c("vazao", "njus")]
    model <- rbindlist(dat)

    if(is.null(vcov)) vcov <- matrix(0, 5, 5)

    out <- list(coefs = coefs, bounds = bounds, model = model, vcov = vcov)
    attr(out, "npoli") <- length(coef)

    fitted <- fitted.polijusU(out)

    out$fitted <- fitted

    attr(out, "tag") <- tag
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

    npoli <- attr(polijus, "npoli")

    coefs <- polijus$coefs

    fitted <- lapply(seq(npoli), function(i) {
        bounds <- polijus$bounds[[i]]
        vaz    <- polijus$model[vazao %between% bounds, vazao]
        coefI  <- coefs[[i]]
        fit    <- lapply(seq(coefI), function(i) vaz^(i - 1) * coefI[i])

        rowSums(do.call(cbind, fit))
    })

    fitted <- unlist(fitted)

    return(fitted)
}

#' @export 

predict.polijusU <- function(polijus, newdata, ...) {

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
        polijus$model[, 1:2 := mapply(.SD, mu, sig, SIMPLIFY = FALSE, FUN = function(d, u, s) d * s + u)]
        polijus$vcov   <- do.call(Matrix::bdiag, DELTA) %*% polijus$vcov %*% do.call(Matrix::bdiag, lapply(DELTA, t))
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
        polijus$model[, 1:2 := mapply(.SD, mu, sig, SIMPLIFY = FALSE, FUN = function(d, u, s) (d - u) / s)]
        polijus$vcov   <- do.call(Matrix::bdiag, DELTA) %*% polijus$vcov %*% do.call(Matrix::bdiag, lapply(DELTA, t))
        polijus$fitted <- fitted(polijus)
    }

    return(polijus)
}