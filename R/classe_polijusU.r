
#' @export

new_polijusU <- function(coefs, bounds, dat, vcov, tag) {

    if(!is.list(coefs)) {
        coefs <- list(poli1 = coefs)
    } else {
        names(coefs) <- paste0("poli", seq(length(coefs)))
    }

    if(!is.list(bounds)) {
        bounds <- list(poli1 = bounds)
    } else {
        names(bounds) <- paste0("poli", seq(length(bounds)))
    }

    dat$hist <- dat$hist[, c("vazao", "njus")]
    model <- rbindlist(dat)

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
        coef <- c(coef, double(5 - length(coef)))

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
    fitted <- fitted[!duplicated(fitted)] # caso haja pontos exatamente iguais aos bounds (pouco provavel)

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