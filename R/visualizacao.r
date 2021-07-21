############################## FUNCOES PARA VISUALIZACAO DE RESULTADOS #############################

#' Plot de objetos \code{datpoli}
#' 
#' Metodo S3 para facilitação dos plots de objetos \code{datpoli}
#' 
#' Exite uma ampla gama de dados contidos em cada objeto \code{datpoli}, de modo que um único 
#' gráfico não é suficiente. Através do parâmetro \code{qual} é possível especificar qual parte do
#' dado deve ser plotada.
#' 
#' Inicialmente, se são desejados plots do dado como um todo, antes de sua classificação por 
#' patamares, é possível fornecer uma das seguintes opções.
#' 
#' \describe{
#' \item{\code{qual} = "todos"}{plota dados brutos, estáveis e filtrados}
#' \item{\code{qual} = "brutos"}{plota apenas dados brutos}
#' \item{\code{qual} = "estaveis"}{plota apenas dados estáveis}
#' \item{\code{qual} = "filtrados"}{plota apenas dados filtrados}
#' }
#' 
#' Deve ser notado que a opção \code{"filtrados"} só é viável após a classificação e filtragem de 
#' dados por patamar através de \code{\link{classfiltrapats}}. Similarmente, a opção \code{"todos"} 
#' só plota os dados filtrados após execução desta função.
#' 
#' A segunda alternativa é o plot de dados por patamar. Neste caso o parâmetro \code{qual} deve ser 
#' fornecido no formato \code{"pat_XXX.X"}. Nesta situação serão plotados apenas os dados referentes
#' ao patamar \code{XXX.X} colorido para indicar aqueles que foram removido ao longo da filtragem.
#' 
#' @param dat objeto \code{datpoli}
#' @param qual string indicando o que deve ser plotado. Ver Detalhes
#' 
#' @return plot das informações requisitadas em \code{qual}
#' 
#' @examples 
#' 
#' # detecta vazoes estaveis
#' dat <- filtravazest(dummydata)
#' 
#' # classifica e filtra por patamares
#' dat <- classfiltrapats(dat)
#' 
#' # plot
#' plot(dat, "todos") # dado completo
#' plot(dat, "estaveis") # apenas aqueles em condicao estavel
#' plot(dat, "pat_020.5") # plot de um patamar especifico
#' 
#' \dontrun{
#' 
#' # tentar plots de informacoes ainda nao geradas resulta em erro
#' plot(dummydata, "filtrados") # dummydata so contem dados brutos, nao existe ainda o historico filtrado
#' }
#' 
#' @export
#' 
#' @rdname importadados

plot.datpoli <- function(dat, qual, ...) {

    if(missing(qual)) qual <- "brutos+estaveis+filtrados"
    if(qual == "todos") qual <- "brutos+estaveis+filtrados"

    if(grepl("(brutos)|(estaveis)|(filtrados)", qual)) plot_func <- plota_datfull
    if(grepl("^pat_", qual)) plot_func <- plota_patfiltro
    if(grepl("^conv_", qual)) plot_func <- plota_patconv
    if(grepl("^remanso$", qual)) plot_func <- plota_remanso

    plot_func(dat, qual)
}

#' @export

plot.datcurvabase <- function(dbase, datorig, ...) {

    dplot1 <- dbase[[1]]
    if(!missing("datorig")) {
        ranges1 <- list(range(datorig[[1]]$vazao, na.rm = TRUE), range(datorig[[1]]$njus, na.rm = TRUE))

        aux <- datorig[[2]][!(datahora %in% dplot1$datahora)]
        aux[, base := FALSE]
    } else {
        ranges1 <- list(range(dbase[[1]]$vazao, na.rm = TRUE), range(dbase[[1]]$njus, na.rm = TRUE))
        aux <- NULL
    }
    dplot1 <- rbind(aux, dplot1)
    cores1 <- ifelse(dplot1$base, "green4", "deepskyblue2")

    dplot2   <- dbase[[2]]
    nquadros <- ceiling(length(dplot2) / 2) * 2
    ranges2  <- list(c(min(dplot1$vazao), max(dplot2[[1]]$vazao)),
                     c(min(dplot1$njus), max(sapply(dplot2, function(d) max(d$njus)))))

    quadros <- matrix(c(rep(1, nquadros), seq(nquadros) + 1), ncol = 4)

    dev.new(width = 30, height = 18)
    layout(quadros)
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.7)
    par(mar = c(5.1, 6.1, 4.1, 2.1))
    plot(dplot1$vazao, dplot1$njus, panel.first = grid(col = "grey85"), col = cores1, pch = 16,
        xlim = ranges1[[1]], ylim = ranges1[[2]],
        xlab = expression("Vazão [m"^3*"/s]"), ylab = "Nível de jusante [m]")
    legend("bottomright", inset = .02, title = "Dados ", cex = 1.5,
        legend = c("históricos", "curva base"), pch = 16,
        col = c("deepskyblue2", "green4"))

    par(mar = c(5.1, 1, 4.1, 1))
    for(ext in names(dplot2)) {
        plot(dplot1[base == TRUE]$vazao, dplot1[base == TRUE]$njus, panel.first = grid(col = "grey85"),
            col = "green4", pch = 16,
            xlim = ranges2[[1]], ylim = ranges2[[2]],
            xlab = "", ylab = "", xaxt = "n", yaxt = "n",
            main = paste0("Extensao via: ", ext))
        lines(dplot2[[ext]][, c("vazao", "njus")], lwd = 3, lty = 2, col = "green4")
    }
    legend("bottomright", inset = .02, cex = 1.5,
        legend = c("Curva base", "extensão"), pch = c(16, NA), lty = c(NA, 2), lwd = c(NA, 3),
        col = "green4")
}

# HELPERS ------------------------------------------------------------------------------------------

plota_datfull <- function(dat, qual) {

    ranges <- list(range(dat[[1]]$vazao, na.rm = TRUE), range(dat[[1]]$njus, na.rm = TRUE))

    qual <- strsplit(qual, "\\+")[[1]]

    # tira estavel ou filtrado de qual se o dado ainda nao passou por essas fases
    if(!attr(dat, "filtravazest") & ("estaveis" %in% qual)) {
        qual <- qual[!grepl("estaveis", qual)]
        warning(paste0("'qual' inclui dados de vazao estaveis, porem esta avaliacao ainda nao foi realizada",
            "\n Use polijus::filtravazest"))
    }
    if(!attr(dat, "classfiltrapats") & ("filtrados" %in% qual)) {
        qual <- qual[!grepl("filtrados", qual)]
        warning(paste0("'qual' inclui dados filtrados, porem esta avaliacao ainda nao foi realizada",
            "\n Use polijus::classfiltrapats"))
    }
    if(length(qual) == 0) stop("Nao ha dados para plotar")

    dplot <- list(copy(dat[[1]]), copy(dat[[2]]))
    colunas <- c("vazao", "njus", "valido", "tipo")

    dplot[[1]][, tipo := "brutos"]
    dplot[[1]] <- dplot[[1]][, ..colunas]

    if(("filtrados" %in% qual)) {
        dplot[[2]] <- dplot[[2]][, .SD, .SDcols = c("vazao", "njus", "valido")]
        dplot[[2]][valido == FALSE, tipo := "estaveis"]
        dplot[[2]][valido == TRUE, tipo := "filtrados"]
    } else if("estaveis" %in% qual) {
        dplot[[2]] <- dplot[[2]][, .SD, .SDcols = c("vazao", "njus", "valido")]
        dplot[[2]][, tipo := "estaveis"]
    } else {
        dplot[[2]] <- NULL
    }

    dplot <- rbindlist(dplot)
    dplot[, tipo := factor(tipo, levels = c("brutos", "estaveis", "filtrados"), ordered = TRUE)]
    setorder(dplot, "tipo")
    dplot <- dplot[tipo %in% qual]

    escala <- c(brutos = "gray90", estaveis = "skyblue1", filtrados = "deepskyblue3")
    cores <- unname(dplot[, escala[tipo]])

    plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = cores, pch = 16,
        xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vazão [m"^3*"/s]"), ylab = "Nível de jusante [m]")
    legend("bottomright", inset = .02, title = "Dados ",
        legend = sub("estaveis", "estáveis", qual), pch = 16,
        col = escala[qual])
}

plota_remanso <- function(dat, ...) {

    if(!attr(dat, "evalremanso")) {
        stop(paste0("'qual' se refere ao efeito de remanso, porem esta analise ainda",
            " nao foi realizada", "\n Use polijus::evalremanso"))
    }

    ranges   <- list(range(dat[[1]]$vazao, na.rm = TRUE), range(dat[[1]]$njus, na.rm = TRUE))

    remanso <- sapply(dat$patinfo, "[[", "remanso")
    vazconv <- sapply(dat$patinfo, "[[", "vazconv")

    dplot <- copy(dat$hist_est)

    # Dados sem efeito de remanso
    dplot[, temremanso := TRUE]
    dplot[pat %in% names(remanso[!remanso]), temremanso := FALSE]
    for(p in names(remanso[remanso])) {
        dplot[(pat == p) & (vazao > vazconv[p]), temremanso := FALSE]
    }
    setorderv(dplot, "temremanso", -1)

    dbar <- dplot[, .("com" = sum(temremanso), "sem" = sum(!temremanso)), by = pat]
    setorder(dbar, pat)
    dbar <- t(data.matrix(dbar))[-1, ]

    layout(matrix(c(1, 2), 2, 1))
    plot(dplot$vazao, dplot$njus, pnael.first = grid(col = "grey85"),
        col = ifelse(dplot$temremanso, "deepskyblue2", "green4"), pch = 16,
        xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vazão [m"^3*"/s]"), ylab = "Nível de jusante [m]",
        main = "Dispersão do efeito de remanso no histórico equivalente filtrado")

    barplot(dbar, , xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", border = NA)
    grid(nx = NA, ny = NULL)
    par(new = TRUE)
    barplot(dbar, col = c("deepskyblue2", "green4"), names.arg = sort(unique(dplot$pat)),
        xlab = "Patamar", ylab = "Contagem", main  = "Número de registros por patamar")

    legend("topleft", inset = 0.02, legend = c("Com remanso", "Sem remanso"),
        fill = c("deepskyblue2", "green4"))
}

plota_patfiltro <- function(dat, qual) {

    if(!attr(dat, "classfiltrapats")) {
        stop(paste0("'qual' se refere aos dados de um patamar especifico, porem a classificacao ainda",
            " nao foi realizada", "\n Use polijus::classfiltrapats"))
    }

    ranges <- list(range(dat[[1]]$vazao, na.rm = TRUE), range(dat[[1]]$njus, na.rm = TRUE))

    qual <- sub("pat_", "", qual)

    dplot  <- copy(dat[[2]])[pat == qual]
    filtro <- dat$patinfo[[qual]]

    vx <- seq(min(dplot$vazao), max(dplot$vazao), length.out = 500)
    ly <- lapply(filtro[[1]], predict, newdata = data.frame(vazao = vx))

    escala <- c("deepskyblue2", "purple", "red", "orange", "yellow2")
    dplot[, filtro := filtro$filtro]
    dplot[, cores := escala[filtro+ 1]]

    plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = dplot$cores, pch = 16,
        xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vazão [m"^3*"/s]"), ylab = "Nível de jusante [m]")
}

plota_patconv <- function(dat, qual) {

    if(!attr(dat, "evalremanso")) {
        stop(paste0("'qual' se refere a convergencia entre dois patamares, porem esta analise ainda",
            " nao foi realizada", "\n Use polijus::evalremanso"))
    }

    ranges <- list(range(dat[[1]]$vazao, na.rm = TRUE), range(dat[[1]]$njus, na.rm = TRUE))

    patamares <- sort(unique(dat$hist_est$pat))

    qual    <- sub("conv_", "", qual)
    qualant <- patamares[which(patamares == qual) - attr(dat, "step_converg")]
    quais   <- c(qualant, qual)

    dplot <- copy(dat[[2]])[pat %in% quais]
    tends <- lapply(dat$patinfo[quais], function(l) l$tend[[4]])

    tends <- mapply(quais, tends, SIMPLIFY = FALSE, FUN = function(pat, mod) {
        vaz  <- if(class(mod) == "loess") mod$x[, 1] else mod$model[, 2]
        njus <- if(class(mod) == "loess") mod$fitted else mod$fitted.values
        out <- data.frame(vazao = vaz, njus = njus, pat = pat)
        out[order(out$vazao), ]
    })

    vazconv     <- dat$patinfo[[qual]]["vazconv"]
    vazconv_reg <- dat$patinfo[[qual]]["vazconv_reg"]

    cores <- structure(c("deepskyblue", "deepskyblue4"), names = quais)

    plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = cores, pch = 16,
        xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vazão [m"^3*"/s]"), ylab = "Nível de jusante [m]")
    for(i in seq(tends)) lines(tends[[i]]$vazao, tends[[i]]$njus, col = cores[i], lwd = 2)
    abline(v = vazconv, lwd = 2, lty = 2, col = "grey85")
    abline(v = vazconv_reg, lwd = 2, lty = 2, col = "black")
    legend("bottomright", inset = .02,
        lwd = 2, col = cores, legend = paste0("Patamar ", quais))
}