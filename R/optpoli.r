################################## OTIMIZACAO DOS HIPERPARAMETROS ##################################

#' Otimizacao de curvas polinomiais por partes
#' 
#' Genérica e métodos para otimizacao dos ajustes polinomiais por partes a cada conjunto de dados
#' 
#' Esta funcao opera de maneria extremamente similar à \code{\link{fitpoli}}, sendo apenas um 
#' wrapper para otimização dos hiperparâmetros do ajuste (i.e. \code{graus}, \code{pto_turbmax} e 
#' \code{pto_ext} ou \code{vaz_ext}).
#' 
#' \code{graus} aqui deve ser fornecido como uma lista de vetores, onde cada elemento da lista
#' contém quais graus devem ser testados para cada parte polinomial (ver Exemplos para mais 
#' detalhes). Como a decisão da configuração de graus é inteira, o que se faz é otimizar os 
#' hiperparâmetros restantes para cada possível combinação de graus fornecida e selecionar aquela
#' cujo resultado final leval ao menor erro de ajuste.
#' 
#' \code{pto_turbmax0}, \code{pto_ext0} ou \code{vaz_ext0} são pontos os valores iniciais para 
#' otimização e devem necessariamente ser supridos.
#' 
#' Por fim, \code{opcoes} pode ser passada como uma lista de até tres elementos lógicos nomeados
#' 
#' \describe{
#'     \item{min_descont}{indicando se a descontinuidade entre partes deve ser minimizada}
#'     \item{fix_turbmax}{indicando se o valor \code{pto_turbmax0} deve ser fixado na otimizacao}
#'     \item{fix_ext}{indicando se o valor \code{(pto)|(vaz)_ext0} deve ser fixado na otimizacao}
#' }
#' 
#' A opção \code{min_descont} só é utilizada em casos de ajuste de dois polinômios aos dados 
#' históricos, onde a suavidade na conexão não é obrigatória. Caso \code{opcoes$min_descont = TRUE}, 
#' após a otimização dos outros hiperparâmetros visando reduzir a descontinuidade da primeira 
#' derivada da curva.
#' 
#' \code{fix_turbmax} e \code{fix_ext} permitem travar os valores iniciais \code{pto_turbmax0},
#' \code{pto_ext0} ou \code{vaz_ext0}
#' 
#' @param dat objeto do tipo \code{datcbase} ou \code{data.table}. Ver Detalhes
#' @param ext dado para extensão do ajuste. Ver Detalhes
#' @param graus lista de até três posições indicando os graus de cada parte polinomial. Ver Detalhes
#' @param pto_turbmax0 vetor de duas posicoes indicando coordenadas da conexão entre polinômios 
#'     ajustados aos dados históricos
#' @param ... demais parâmetros específicos de cada método. Ver Detalhes
#' @param opcoes lista contendo opcoes extras para a otimizacao. Ver Detalhes
#' 
#' @return lista de objetos \code{polijusU}. Ver detalhes
#' 
#' @seealso \code{\link{polijusU}} para detalhes do objeto retornado; \code{\link{fitpoli}} para 
#'     execucao de ajuste simples
#' 
#' @export

optpoli <- function(dat, ext, graus, pto_turbmax0, ..., opcoes) UseMethod("optpoli")

#' @param pto_ext0 vetor de duas posicoes indicando coordenadas da conexão entre último polinômio
#'     ajustados aos dados históricos e polinômio dos dados de extensão
#' 
#' @examples 
#' 
#' # extrai uma curva base
#' dbase <- extraibase(dummydata)
#' 
#' # selecao do dado de extensao
#' ext <- 1
#' 
#' # um polinomio ao historico e um para extensao ---------------------------------------
#' 
#' # a lista de duas posicoes indica o ajuste de dois polinomios, testando os graus 2, 3 e 4
#' # em cada parte, para um total de 9 configuracoes possiveis
#' graus <- list(2:4, 2:4)
#' 
#' # heuristica para estimar a inicializacao do ponto de conexao com o dado de extensao
#' maxvaz   <- 1.1 * dbase$hist[, max(vazao)]
#' pto_ext0 <- c(maxvaz, dummydata$ext[[ext]](maxvaz))
#' 
#' ajustes2poli <- optpoli(dbase, "CAD", graus = graus, pto_ext0 = pto_ext0)
#' 
#' # dois polinomios historicos e um de extensao ----------------------------------------
#' 
#' # aumentamos a lista de graus para testar, pois agora sao tres polinomios
#' graus <- list(2:4, 2:4, 2:4)
#' 
#' # heuristica para estimar a inicializacao do ponto de conexao entre polinomios do dado historico
#' engef    <- attr(dbase, "vazef")
#' yturbmax <- dbase$hist[(vazao < engef * 1.05) & (vazao > engef * 0.95), mean(njus)]
#' pto_turbmax0 <- c(engef, yturbmax)
#' 
#' # utilizando minimizacao de descontinuidade
#' ajustes3poli <- optpoli(dbase, "CAD", graus = graus, pto_ext0 = pto_ext0,
#'     pto_turbmax0 = pto_turbmax0, opcoes = list(min_descont = TRUE))
#' 
#' # por fim, podemos ver uma compilacao de todos os ajustes com
#' do.call(summary, ajustes3poli)
#' 
#' @rdname optpoli
#' 
#' @export

optpoli.datcbase <- function(dat, ext, graus, pto_turbmax0, pto_ext0, opcoes, ...) {

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

                optpar <- optim(optpar, obj2, control = list(maxit = 200, alpha = -1))$par

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

    erros <- sapply(out, function(x) mean(unlist(residuals(x))^2))

    out <- out[order(erros)]

    return(out)
}

#' @param vaz_ext0 escalar indicando a vazão da conexão entre último polinômio ajustados aos dados
#'     históricos e curva base
#' @param zero_forcado booleano indicando se o A0 deve ser restrito igual ao patamar de referência
#' 
#' @examples 
#' 
#' # otimizacao de ajuste individual ----------------------------------------------------
#' 
#' # utilizando a segunda melhor curva base do ajuste2poli
#' polibase <- ajustes2poli[[2]]
#' 
#' # utilizando patamer 33.0
#' datind <- extraipats(dummydata, quais = 33)[[1]]
#' 
#' # inicializacao de vaz_ext0
#' vaz_ext0 <- datind[, max(vazao)]
#' 
#' # otimizacao com apenas um polinomio
#' ajustesind <- optpoli(datind, polibase, list(2:4), vaz_ext0 = vaz_ext0)
#' 
#' do.call(summary, ajustesind)
#' 
#' @rdname optpoli
#' 
#' @export

optpoli.default <- function(dat, ext, graus, pto_turbmax0, vaz_ext0, zero_forcado, opcoes, ...) {

    vazao <- njus <- NULL

    graus <- data.matrix(do.call(expand.grid, graus))

    if(missing("ext")) ext <- c(NA)
    if(missing("pto_turbmax0")) pto_turbmax0 <- c(NA, NA)
    if(missing("vaz_ext0")) vaz_ext0 <- NA
    if(missing("zero_forcado")) zero_forcado <- NA

    if(missing("opcoes")) {
        opcoes <- list(min_descont = FALSE, fix_turbmax = FALSE, fix_ext = FALSE)
    } else {
        if(is.null(opcoes$min_descont)) opcoes$min_descont <- FALSE
        if(is.null(opcoes$fix_turbmax)) opcoes$fix_turbmax <- FALSE
        if(is.null(opcoes$fix_ext)) opcoes$fix_ext <- FALSE
    }

    opcoes$min_descont <- opcoes$min_descont & !all(is.na(pto_turbmax0))
    opcoes$fix_turbmax <- opcoes$fix_turbmax | all(is.na(pto_turbmax0))
    opcoes$fix_ext     <- opcoes$fix_ext | all(is.na(vaz_ext0))

    if(!(opcoes$fix_turbmax | opcoes$fix_ext)) {
        mapfunc <- function(x, func) {func$pto_turbmax <- x[1:2]; func$vaz_ext <- x[3]; func}
    } else if(opcoes$fix_turbmax & opcoes$fix_ext) {
        mapfunc <- function(x, func) func
    } else if(opcoes$fix_ext) {
        mapfunc <- function(x, func) {func$pto_turbmax <- x[1:2]; func}
    } else if(opcoes$fix_turbmax) {
        mapfunc <- function(x, func) {func$vaz_ext <- x[1]; func}
    }

    cli::cli_h2("Otimizacao de hiperparametros")

    out <- lapply(seq(nrow(graus)), function(i) {
        l_parse <- parseargsind(dat, ext, graus[i, ], pto_turbmax0, vaz_ext0, zero_forcado)
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

            par0 <- unlist(l_func[c("pto_turbmax", "vaz_ext")])
            optpar <- optim(par0[!is.na(par0)], obj1, method = "BFGS", control = list(maxit = 200))$par

            l_func <- mapfunc(optpar, l_func)

            fit <- eval(do.call(call, l_func))
            fit <- rescale(fit, scales)

            if(opcoes$min_descont) {

                # ver comentario em optpoli.datcbase

                optpar <- optim(optpar, obj2, control = list(maxit = 200, alpha = -1))$par

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

    erros <- sapply(out, function(x) mean(unlist(residuals(x))^2))

    out <- out[order(erros)]

    return(out)
}