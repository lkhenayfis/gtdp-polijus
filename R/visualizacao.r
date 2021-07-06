############################## FUNCOES PARA VISUALIZACAO DE RESULTADOS #############################

#' @import data.table
#' 
#' @export

plot.datpoli <- function(dat, tipo = c("bruto", "estavel", "filtrado"), ...) {

    dplot <- list(copy(dat[[1]]), copy(dat[[2]]))
    colunas <- c("vazao", "njus", "valido", "tipo")

    dplot[[1]][, tipo := "brutos"]
    dplot[[1]] <- dplot[[1]][, ..colunas]

    if(attr(dat, "estavel") == TRUE) {
        if(attr(dat, "filtrado") == TRUE) {
            dplot[[2]] <- dplot[[2]][, .SD, .SDcols = c("vazao", "njus", "valido")]
            dplot[[2]][valido == FALSE, tipo := "estaveis"]
            dplot[[2]][valido == TRUE, tipo := "filtrados"]
        } else {
            dplot[[2]] <- dplot[[2]][, .SD, .SDcols = c("vazao", "njus", "valido")]
            dplot[[2]][, tipo := "estaveis"]
        }
    }

    dplot <- rbindlist(dplot)
    dplot[, tipo := factor(tipo, levels = unique(tipo))]

    escala <- c(brutos = "slategray2", estaveis = "deepskyblue2", filtrados = "skyblue2")
    cores <- unname(dplot[, escala[tipo]])

    plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = cores,
        xlab = "Vazão [m³/s]", ylab = "Nível de jusante [m]")
    legend("bottomright", inset = .02, title = "Dados ",
        pch = 1, legend = levels(dplot$tipo),
        col = escala)
}
