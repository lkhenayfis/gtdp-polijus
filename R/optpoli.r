################################## OTIMIZACAO DOS HIPERPARAMETROS ##################################

optpoli <- function(dat, ext, graus, pto_turbmax0, ...) UseMethod("optpoli")

optpoli.datcbase <- function(dat, ext, graus, pto_turbmax0, pto_ext0, opcoes) {

    vazao <- njus <- NULL

    graus <- data.matrix(do.call(expand.grid, graus))

    if(missing("ext")) ext <- c(NA)
    if(missing("pto_turbmax0")) pto_turbmax0 <- c(NA, NA)
    if(missing("pto_ext0")) pto_ext0 <- c(NA, NA)

    if(missing("opcoes")) {
        opcoes <- list(min_descont = FALSE, fix_turbmax = FALSE, fix_ext = FALSE)
    } else {
        if(is.null(opcoes$min_descont)) opcoes$min_descont <- FALSE
        if(is.null(opcoes$fix_turbmax)) opcoes$fix_turbmax <- FALSE
        if(is.null(opcoes$fix_ext)) opcoes$fix_ext <- FALSE
    }

    opcoes$min_descont <- opcoes$min_descont & !all(is.na(pto_turbmax0))
    opcoes$fix_turbmax <- opcoes$fix_turbmax | all(is.na(pto_turbmax0))
    opcoes$fix_ext     <- opcoes$fix_ext | all(is.na(pto_ext0))

    if(!(opcoes$fix_turbmax | opcoes$fix_ext)) {
        mapfunc <- function(x, func) {func$pto_turbmax <- x[1:2]; func$pto_ext <- x[3:4]; func}
    } else if(opcoes$fix_turbmax & opcoes$fix_ext) {
        mapfunc <- function(x, func) func
    } else if(opcoes$fix_ext) {
        mapfunc <- function(x, func) {func$pto_turbmax <- x[1:2]; func}
    } else if(opcoes$fix_turbmax) {
        mapfunc <- function(x, func) {func$pto_ext <- x[1:2]; func}
    }

    cli::cli_h2("Otimizacao de hiperparametros")

    out <- lapply(seq(nrow(graus)), function(i) {
        l_parse <- parseargsbase(dat, ext, graus[i, ], pto_turbmax0, pto_ext0)
        l_func  <- l_parse[[1]]
        scales  <- l_parse[[2]]

        obj1 <- function(x) {
            l_func <- mapfunc(x, l_func)
            fit <- eval(do.call(call, l_func))
            fit <- rescale(fit, scales)

            mean(unlist(residuals(fit))^2)
        }

        obj2 <- function(x) {
            l_func <- mapfunc(x, l_func)
            fit <- eval(do.call(call, l_func))
            fit <- rescale(fit, scales)

            new <- data.frame(vazao = c(x[1] - 1e-5, x[1] + 1e-5))
            diff(predict(fit, newdata = new, derivada = TRUE))^2
        }

        tryCatch({

            par0 <- unlist(l_func[c("pto_turbmax", "pto_ext")])
            optpar <- optim(par0[!is.na(par0)], obj1, control = list(maxit = 2))$par

            l_func <- mapfunc(optpar, l_func)

            fit <- eval(do.call(call, l_func))
            fit <- rescale(fit, scales)

            if(opcoes$min_descont) {

                # duas observacoes devem ser feitas aqui. Primeiro, tanto o pronto de conexao entre
                # polinomios historicos quanto entre o segundo e extensao sao modeificados nessa
                # parte. Isto e feito porque ambos os pontos de controle interferem no melhor ajuste
                # do polinomio
                # Em segundo lugar, alpha = -1 para impedir que o nelder mead faca reflexao (e assim
                # tambem nao faz expansao) ficando so a parte de contracao. O objetivo dessa
                # parametrizacao e evitar que a segunda otimizacao fuja muito do primeiro optpar

                optpar <- optim(optpar, obj2, control = list(maxit = 2, alpha = -1))$par

                l_func <- mapfunc(optpar, l_func)

                fit <- eval(do.call(call, l_func))
                fit <- rescale(fit, scales)
            }

            cli::cli_alert_success(paste0("Graus: (", paste0(graus[i, ], collapse = ", "), ")"))
            return(fit)
        }, error = function(e) {

            cli::cli_alert_danger(paste0("Graus: (", paste0(graus[i, ], collapse = ", "), ")"))
            return(NULL)
        })
    })

    out <- out[!sapply(out, is.null)]

    return(out)
}