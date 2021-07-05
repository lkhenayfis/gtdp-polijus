############################# LEITURA E FORMATACAO DOS DADOS DE ENTRADA ############################

#' @import data.table
#' @importFrom readxl read_xlsx

importadados <- function(path = NULL) {

    if(is.null(path)) path <- file.choose()

    # Le cadastro e pega codigo da usina
    cad <- read_xlsx(path, sheet = 1)
    cod <- as.numeric(cad[1, 2])

    vazef    <- HIDR[HIDR$USI == cod, "VAZEF"]
    policoef <- unlist(HIDR[HIDR$USI == cod, -1])

    # Leitura do historico
    hist <- setDT(read_xlsx(path, sheet = 2))
    colnames(hist) <- c("data", "hora", "njus", "vazao", "nmont", "valido")
    hist <- tratahist(hist)

    # Leitura dos dados de extensao
    ext <- setDT(read_xlsx(path, sheet = "Extensao", col_types = "numeric", skip = 1))
    ext <- lapply(seq(ncol(ext))[c(TRUE, FALSE)], function(i) {
        cols <- c(i, i + 1)
        out <- ext[, ..cols]
        colnames(out) <- c("vazao", "njus")
        out
    })
    ext <- trataext(ext)

}

#' @import data.table

tratahist <- function(hist) {

    if(!is.data.table(hist)) hist <- as.data.table(hist)

    colunas <- c("data", "hora", "njus", "vazao", "nmont", "valido")
    if(!all(colunas %in% colnames(hist))) {
        stop("Faltam as colunas\n", colnames(hist)[!(colunas %in% colnames(hist))])
    } else {
        setcolorder(hist, colunas)
    }

    hist[, datahora := paste(data, paste0((hora - 1), ":00:00"))]
    hist[, datahora := fasttime::fastPOSIXct(datahora, tz = "GMT")]
    if(!all(diff(hist$datahora) == 1)) stop("Ha datas faltantes no historico")

    hist[vazao < 0, valido := FALSE]
    hist[njus < 0, valido := FALSE]
    hist[nmont < 0, valido := FALSE]

    return(hist)
}

#' @import data.table

trataext <- function(ext) {

    if(class(ext) != "list") ext <- list(ext)
    for(i in seq(ext)) {
        ext[[i]] <- as.data.table(ext[[i]])
        if(!all(c("vazao", "njus") %in% colnames(ext[[i]]))) stop("Extensao #", i, " nao contem colunas necessarias")
    }
    if(length(ext) > 0) names(ext) <- paste0("EXT", seq(ext))

    ext <- ext[!sapply(ext, function(d) all(is.na(d)))]

    ext <- lapply(ext, function(d) {
        d   <- setorder(d, "vazao")
        mod <- loess(njus ~ vazao, d)
        function(x) predictCpp2(d$vazao, mod$fitted, x)
    })

    return(ext)
}