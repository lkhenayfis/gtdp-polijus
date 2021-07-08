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

    patamar <- grepl("pat_", qual)

    ranges <- list(range(dat[[1]]$vazao, na.rm = TRUE), range(dat[[1]]$njus, na.rm = TRUE))

    if(!patamar) {

        qual <- strsplit(qual, "\\+")[[1]]

        # tira estavel ou filtrado de qual se o dado ainda nao passou por essas fases
        if(!attr(dat, "estavel") & ("estaveis" %in% qual)) {
            qual <- qual[!grepl("estaveis", qual)]
            warning(paste0("'qual' inclui dados de vazao estaveis, porem esta avaliacao ainda nao foi realizada",
                "\n Use polijus::filtravazest()"))
        }
        if(!attr(dat, "classpats") & ("filtrados" %in% qual)) {
            qual <- qual[!grepl("filtrados", qual)]
            warning(paste0("'qual' inclui dados filtrados, porem esta avaliacao ainda nao foi realizada",
                "\n Use polijus::classfiltrapats()"))
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
        }

        dplot <- rbindlist(dplot)
        dplot[, tipo := factor(tipo, levels = c("brutos", "estaveis", "filtrados"), ordered = TRUE)]
        setorder(dplot, "tipo")
        dplot <- dplot[tipo %in% qual]

        escala <- c(brutos = "gray90", estaveis = "skyblue1", filtrados = "deepskyblue3")
        cores <- unname(dplot[, escala[tipo]])

        plot(dplot$vazao, dplot$njus, panel.first = grid(col = "grey85"), col = cores, pch = 16,
            xlim = ranges[[1]], ylim = ranges[[2]],
            xlab = "Vazão defluente [m³/s]", ylab = "Nível de jusante [m]")
        legend("bottomright", inset = .02, title = "Dados ",
            legend = sub("estaveis", "estáveis", qual), pch = 16,
            col = escala[qual])
    } else {

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
            xlab = "Vazão defluente [m³/s]", ylab = "Nível de jusante [m]")
    }
}
