################################## OTIMIZACAO DOS HIPERPARAMETROS ##################################

optpoli <- function(dat, ext, graus, pto_turbmax0, ...) UseMethod("optpoli")

optpoli.datcbase <- function(dat, ext, graus, pto_turbmax0, pto_ext0, opcoes) {

    vazao <- njus <- NULL

    graus <- data.matrix(do.call(expand.grid, graus))

    if(missing("opcoes")) {
        opcoes <- list(min.descont = FALSE, fix_turbmax = FALSE, fix_ext = FALSE)
    } else {
        if(is.null(opcoes$min_descont)) opcoes$min_descont <- FALSE
        if(is.null(opcoes$fix_turbmax)) opcoes$fix_turbmax <- FALSE
        if(is.null(opcoes$fix_ext)) opcoes$fix_ext <- FALSE
    }

    opcoes$fix_turbmax <- opcoes$fix_turbmax & missing("pto_turbmax0")
    opcoes$fix_ext     <- opcoes$fix_ext & missing("pto_ext0")

    if(!(opcoes$fix_turbmax | opcoes$fix_ext)) {
        par0 <- c(pto_turbmax0, pto_ext0)
        mapfunc <- function(x, func) {func$pto_turbmax <- x[1:2]; func$pto_ext <- x[3:4]; func}
    } else if(opcoes$fix_turbmax & opcoes$fix_ext) {
        par0 <- numeric(0)
        mapfunc <- function(x, func) func
    } else if(opcoes$fix_ext) {
        par0 <- pto_turbmax0
        mapfunc <- function(x, func) {func$pto_turbmax <- x[1:2]; func}
    } else if(opcoes$fix_turbmax) {
        par0 <- pto_ext0
        mapfunc <- function(x, func) {func$pto_ext <- x[1:2]; func}
    }

    cli::cli_h2("Otimizacao de hiperparametros")

    out <- lapply(seq(nrow(graus)), function(i) {
        l_parse <- parseargsbase(dat, ext, graus[i, ], pto_turbmax0, pto_ext0)
        l_func  <- l_parse[[1]]
        scales  <- l_parse[[2]]

        obj <- function(x) {
            l_func <- mapfunc(x, l_func)
            fit <- eval(do.call(call, l_func))
            fit <- rescale(fit, scales)

            mean(unlist(residuals(fit))^2)
        }

        tryCatch({

            par0 <- unlist(l_func[c("pto_turbmax", "pto_ext")])
            optpar <- optim(par0, obj, control = list(maxit = 200))$par

            l_func <- mapfunc(optpar, l_func)

            fit <- eval(do.call(call, l_func))
            fit <- rescale(fit, scales)

            cli::cli_alert_success(paste0("Graus: (", paste0(graus[i, ], collapse = ", "), ")"))
            return(fit)
        }, error = function(e) {

            cli::cli_alert_danger(paste0("Graus: (", paste0(graus[i, ], collapse = ", "), ")"))
            return(NULL)
        })
    })

    return(out)
}