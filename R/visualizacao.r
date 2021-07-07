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
#' @param dat objeto \code{datpolo}
#' @param qual string indicando o que deve ser plotado. Ver Detalhes
#' 
#' @return plot das informações requisitadas em \code{qual}
#' 
#' @export
#' 
#' @rdname importadados

plot.datpoli <- function(dat, qual, ...) {

    if(missing(qual)) qual <- "brutos+estaveis+filtrados"
    if(qual == "todos") qual <- "brutos+estaveis+filtrados"

    patamar <- grepl("pat_", qual)

    if(!patamar) {

        qual <- strsplit(qual, "\\+")[[1]]

        dplot <- list(copy(dat[[1]]), copy(dat[[2]]))
        colunas <- c("vazao", "njus", "valido", "tipo")

        dplot[[1]][, tipo := "brutos"]
        dplot[[1]] <- dplot[[1]][, ..colunas]

        if(attr(dat, "estavel") == TRUE) {
            if(attr(dat, "classpats") == TRUE) {
                dplot[[2]] <- dplot[[2]][, .SD, .SDcols = c("vazao", "njus", "valido")]
                dplot[[2]][valido == FALSE, tipo := "estaveis"]
                dplot[[2]][valido == TRUE, tipo := "filtrados"]
            } else {
                dplot[[2]] <- dplot[[2]][, .SD, .SDcols = c("vazao", "njus", "valido")]
                dplot[[2]][, tipo := "estaveis"]
            }
        }

        dplot <- rbindlist(dplot)
        dplot[, tipo := factor(tipo, levels = c("brutos", "estaveis", "filtrados"), ordered = TRUE)]
        setorder(dplot, "tipo")
        dplot <- dplot[tipo %in% qual]

        escala <- c(brutos = "slategray2", estaveis = "deepskyblue4", filtrados = "purple2")
        cores <- unname(dplot[, escala[tipo]])

        plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = cores,
            xlab = "Vazão defluente [m³/s]", ylab = "Nível de jusante [m]")
        legend("bottomright", inset = .02, title = "Dados ",
            pch = 1, legend = levels(dplot$tipo),
            col = escala)
    } else {

        qual <- sub("pat_", "", qual)

        dplot  <- copy(dat[[2]])[pat == qual]
        filtro <- dat$patinfo[[qual]]

        vx <- seq(min(dplot$vazao), max(dplot$vazao), length.out = 500)
        ly <- lapply(filtro[[1]], predict, newdata = data.frame(vazao = vx))

        ranges <- list(range(dat[[1]]$vazao, na.rm = TRUE), range(dat[[1]]$njus, na.rm = TRUE))

        escala <- c("deepskyblue2", "purple", "red", "orange", "yellow2")
        dplot[, filtro := filtro$filtro]
        dplot[, cores := escala[filtro+ 1]]

        plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = dplot$cores,
            xlim = ranges[[1]], ylim = ranges[[2]],
            xlab = "Vazão defluente [m³/s]", ylab = "Nível de jusante [m]")
    }

}
