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
#' @param x objeto \code{datpoli}
#' @param qual string indicando o que deve ser plotado. Ver Detalhes
#' @param legenda booleano indicando se a legenda deve ou nao ser adicionada
#' @param ... não é utilizado, existe apenas para compatibilização com a genérica \code{plot}
#' 
#' @return plot das informações requisitadas em \code{qual}
#' 
#' @examples 
#' 
#' plot(dummydata, "todos") # dado completo
#' plot(dummydata, "estaveis") # apenas aqueles em condicao estavel
#' plot(dummydata, "pat_020.5") # plot de um patamar especifico
#' 
#' @export

plot.datpoli <- function(x, qual, legenda = TRUE, ...) {

    if(missing(qual)) qual <- "brutos+estaveis+filtrados"
    if(qual == "todos") qual <- "brutos+estaveis+filtrados"

    if(grepl("(brutos)|(estaveis)|(filtrados)", qual)) plot_func <- plota_datfull
    if(grepl("^pat_", qual)) plot_func <- plota_patfiltro
    if(grepl("^conv_", qual)) plot_func <- plota_patconv
    if(grepl("^remanso$", qual)) plot_func <- plota_remanso

    plot_func(x, qual, legenda = legenda)
}

#' @export

plot.datcbase <- function(x, datorig, ...) {

    datahora <- valido <- base <- NULL

    dplot1 <- x[[1]]
    if(!missing("datorig")) {
        ranges1 <- datorig[[1]][valido == TRUE, lapply(.SD, range), .SDcols = c("vazao", "njus")]
        ranges1 <- as.list(ranges1)

        aux <- datorig[[2]][!(datahora %in% dplot1$datahora) & (valido == TRUE)]
        aux[, base := FALSE]
    } else {
        ranges1 <- list(range(x[[1]]$vazao, na.rm = TRUE), range(x[[1]]$njus, na.rm = TRUE))
        aux <- NULL
    }
    dplot1 <- rbind(aux, dplot1)
    cores1 <- ifelse(dplot1$base, "green4", "deepskyblue2")

    dplot2   <- x[[2]]
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
        xlab = expression("Vaz\u00E3o [m"^3 * "/s]"), ylab = "N\u00EDvel de jusante [m]")
    legend("bottomright", inset = .02, title = "Dados ", cex = 1.5,
        legend = c("hist\u00F3ricos", "curva base"), pch = 16,
        col = c("deepskyblue2", "green4"))

    par(mar = c(5.1, 1, 4.1, 1))
    for(ext in names(dplot2)) {
        plot(dplot1[base == TRUE]$vazao, dplot1[base == TRUE]$njus, panel.first = grid(col = "grey85"),
            col = "green4", pch = 16,
            xlim = ranges2[[1]], ylim = ranges2[[2]],
            xlab = "", ylab = "", xaxt = "n", yaxt = "n",
            main = paste0("Extens\u00E3o via: ", ext))
        lines(dplot2[[ext]][, c("vazao", "njus")], lwd = 3, lty = 2, col = "green4")
    }
    legend("bottomright", inset = .02, cex = 1.5,
        legend = c("Curva base", "extens\u00E3o"), pch = c(16, NA), lty = c(NA, 2), lwd = c(NA, 3),
        col = "green4")
}

#' @export

plot.polijusU <- function(x, ...) {

    npoli <- attr(x, "npoli")

    bounds <- unlist(x$bounds)

    vx <- seq(head(bounds, 1), tail(bounds, 1), length.out = 500 * npoli)
    vy <- predict(x, newdata = data.frame(vazao = vx))
    poli <- findInterval(vx, unique(bounds), all.inside = TRUE)

    cores <- c("gold1", "orange1", "coral")

    titulo <- attr(x, "tag")
    if(grepl("curva base", titulo)) {
        titulo <- paste0("Ajuste da ", titulo)
    } else {
        titulo <- paste0("Ajuste do ", titulo)
    }

    rangex <- c(min(min(vx), min(x$model$vazao)), max(max(vx), max(x$model$vazao)))
    rangey <- c(min(min(vy), min(x$model$njus)),  max(max(vy), max(x$model$njus)))

    plot(x$model[, c("vazao", "njus")], panel.first = grid(col = "grey85"),
        col = "deepskyblue2", pch = 16, xlim = rangex, ylim = rangey,
        xlab = expression("Vaz\u00E3o [m"^3 * "/s]"), ylab = "N\u00EDvel de jusante [m]",
        main = titulo)
    for(i in seq(npoli)) {
        lines(vx[poli == i], vy[poli == i], col = cores[i], lwd = 2)
    }
    legend("bottomright", inset = .02,
        legend = c("Dados ajustados", paste0("Polin\u00F4mio ", seq(npoli))),
        pch = c(16, rep(NA, npoli)), lty = c(NA, rep(1, npoli)), lwd = c(NA, rep(2, npoli)),
        col = c("deepskyblue2", cores[seq(npoli)]))
}

#' @export

plot.datpoliind <- function(x, datorig, ...) {

    plot(datorig, "filtrados", FALSE)

    cores <- c("red", "purple", "gold3", "orange3")
    for(i in seq_along(x)) {
        points(x[[i]][, 3:2], col = cores[i], pch = 16)
    }
    legend("bottomright", inset = .02, ncol = 2,
        legend = c("Dados filtrados", paste0("Patamar ", names(x))),
        pch = 16, col = c("deepskyblue1", cores))
}

#' @export

plot.polijusM <- function(x, full = FALSE, ...) {

    ncurvas <- attr(x, "ncurvas")

    cores <- c("green4", "red", "purple", "gold3", "orange3")[seq((ncurvas))]

    if(full) {
        ranges <- list(range(x$curvas[[1]]$model$vazao), range(x$curvas[[1]]$model$njus))
    } else {
        ranges <- list(range(x$model$vazao), range(x$model$njus))
    }
    vx <- x$curvas[[1]]$model$vazao

    plot(x$model$vazao, x$model$njus, panel.first = grid(col = "grey85"), col = "deepskyblue1",
        pch = 16, xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vaz\u00E3o [m"^3 * "/s]"), ylab = "N\u00EDvel de jusante [m]")
    for(i in rev(seq(ncurvas))) {
        lines(vx, predict(x$curvas[[i]], newdata = data.frame(vazao = vx)), col = cores[i], lwd = 2)
    }
    legend("bottomright", inset = .02, ncol = 2,
        legend = c("Dados filtrados", sapply(x$curvas, attr, "tag")),
        pch = c(16, rep(NA, ncurvas)), lwd = c(NA, rep(2, ncurvas)), lty = c(NA, rep(1, ncurvas)),
        col = c("deepskyblue1", cores))
}

# HELPERS ------------------------------------------------------------------------------------------

plota_datfull <- function(x, qual, legenda) {

    valido <- tipo <- NULL

    ranges <- x[[1]][valido == TRUE, lapply(.SD, range), .SDcols = c("vazao", "njus")]
    ranges <- as.list(ranges)

    qual <- strsplit(qual, "\\+")[[1]]

    # tira estavel ou filtrado de qual se o dado ainda nao passou por essas fases
    if(!attr(x, "filtravazest") & ("estaveis" %in% qual)) {
        qual <- qual[!grepl("estaveis", qual)]
        warning(paste0("'qual' inclui dados de vazao estaveis, porem esta avaliacao ainda nao foi realizada",
            "\n Use polijus::filtravazest"))
    }
    if(!attr(x, "classfiltrapats") & ("filtrados" %in% qual)) {
        qual <- qual[!grepl("filtrados", qual)]
        warning(paste0("'qual' inclui dados filtrados, porem esta avaliacao ainda nao foi realizada",
            "\n Use polijus::classfiltrapats"))
    }
    if(length(qual) == 0) stop("Nao ha dados para plotar")

    dplot <- list(copy(x[[1]]), copy(x[[2]]))
    colunas <- c("vazao", "njus", "valido", "tipo")

    dplot[[1]][, tipo := "brutos"]
    dplot[[1]] <- dplot[[1]][, .SD, .SDcols = colunas]

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
        xlab = expression("Vaz\u00E3o [m"^3 * "/s]"), ylab = "N\u00EDvel de jusante [m]")
    if(legenda) {
        legend("bottomright", inset = .02, title = "Dados ",
            legend = sub("estaveis", "est\u00E1veis", qual), pch = 16,
            col = escala[qual])
    }
}

plota_remanso <- function(x, ...) {

    vazao <- pat <- valido <- temremanso <- NULL

    if(!attr(x, "evalremanso")) {
        stop(paste0("'qual' se refere ao efeito de remanso, porem esta analise ainda",
            " nao foi realizada", "\n Use polijus::evalremanso"))
    }

    ranges <- x[[1]][valido == TRUE, lapply(.SD, range), .SDcols = c("vazao", "njus")]
    ranges <- as.list(ranges)

    remanso <- sapply(x$patinfo, "[[", "remanso")
    vazconv <- sapply(x$patinfo, "[[", "vazconv")

    dplot <- copy(x$hist_est)

    # Dados sem efeito de remanso
    dplot[, temremanso := TRUE]
    dplot[pat %in% names(remanso[!remanso]), temremanso := FALSE]
    for(p in names(remanso[remanso])) {
        dplot[(pat == p) & (vazao > vazconv[p]), temremanso := FALSE]
    }
    setorderv(dplot, "temremanso", -1)

    dbar <- dplot[, list("com" = sum(temremanso), "sem" = sum(!temremanso)), by = pat]
    setorder(dbar, pat)
    dbar <- t(data.matrix(dbar))[-1, ]

    layout(matrix(c(1, 2), 2, 1))
    plot(dplot$vazao, dplot$njus, pnael.first = grid(col = "grey85"),
        col = ifelse(dplot$temremanso, "deepskyblue2", "green4"), pch = 16,
        xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vaz\u00E3o [m"^3 * "/s]"), ylab = "N\u00EDvel de jusante [m]",
        main = "Dispers\u00E3o do efeito de remanso no hist\u00F3rico equivalente filtrado")

    barplot(dbar, , xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", border = NA)
    grid(nx = NA, ny = NULL)
    par(new = TRUE)
    barplot(dbar, col = c("deepskyblue2", "green4"), names.arg = sort(unique(dplot$pat)),
        xlab = "Patamar", ylab = "Contagem", main  = "N\u00FAmero de registros por patamar")

    legend("topleft", inset = 0.02, legend = c("Com remanso", "Sem remanso"),
        fill = c("deepskyblue2", "green4"))
}

plota_patfiltro <- function(x, qual, ...) {

    pat <- valido <- filtro <- cores <- NULL

    if(!attr(x, "classfiltrapats")) {
        stop(paste0("'qual' se refere aos dados de um patamar especifico, porem a classificacao ainda",
            " nao foi realizada", "\n Use polijus::classfiltrapats"))
    }

    ranges <- x[[1]][valido == TRUE, lapply(.SD, range), .SDcols = c("vazao", "njus")]
    ranges <- as.list(ranges)

    qual <- sub("pat_", "", qual)

    dplot  <- copy(x[[2]])[pat == qual]
    filtro <- x$patinfo[[qual]]

    escala <- c("deepskyblue2", "purple", "red", "orange", "yellow2")
    dplot[, filtro := filtro$filtro]
    dplot[, cores := escala[filtro + 1]]

    plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = dplot$cores, pch = 16,
        xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vaz\u00E3o [m"^3 * "/s]"), ylab = "N\u00EDvel de jusante [m]")
}

plota_patconv <- function(x, qual) {

    pat <- valido <- NULL

    if(!attr(x, "evalremanso")) {
        stop(paste0("'qual' se refere a convergencia entre dois patamares, porem esta analise ainda",
            " nao foi realizada", "\n Use polijus::evalremanso"))
    }

    ranges <- x[[1]][valido == TRUE, lapply(.SD, range), .SDcols = c("vazao", "njus")]
    ranges <- as.list(ranges)

    patamares <- sort(unique(x$hist_est$pat))

    qual    <- sub("conv_", "", qual)
    qualant <- patamares[which(patamares == qual) - attr(x, "step_converg")]
    quais   <- c(qualant, qual)

    dplot <- copy(x[[2]])[pat %in% quais]
    tends <- lapply(x$patinfo[quais], function(l) l$tend[[4]])

    tends <- mapply(quais, tends, SIMPLIFY = FALSE, FUN = function(pat, mod) {
        vaz  <- if(class(mod) == "loess") mod$x[, 1] else mod$model[, 2]
        njus <- if(class(mod) == "loess") mod$fitted else mod$fitted.values
        out <- data.frame(vazao = vaz, njus = njus, pat = pat)
        out[order(out$vazao), ]
    })

    vazconv     <- x$patinfo[[qual]]["vazconv"]
    vazconv_reg <- x$patinfo[[qual]]["vazconv_reg"]

    cores <- structure(c("deepskyblue", "deepskyblue4"), names = quais)

    plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = cores, pch = 16,
        xlim = ranges[[1]], ylim = ranges[[2]],
        xlab = expression("Vaz\u00E3o [m"^3 * "/s]"), ylab = "N\u00EDvel de jusante [m]")
    for(i in seq(tends)) lines(tends[[i]]$vazao, tends[[i]]$njus, col = cores[i], lwd = 2)
    abline(v = vazconv, lwd = 2, lty = 2, col = "grey85")
    abline(v = vazconv_reg, lwd = 2, lty = 2, col = "black")
    legend("bottomright", inset = .02,
        lwd = 2, col = cores, legend = paste0("Patamar ", quais))
}
