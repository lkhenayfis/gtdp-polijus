path <- "X:/Ciclo 1 - 2005 a 2014/_Resultados Finais por Usina"
dirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)
dirs <- paste0(dirs, "/Niveis de Jusante")
arqs <- file.path(dirs, "Resultados/img.RData")

temresultado <- file.exists(arqs)
arqs <- arqs[temresultado]
arqs <- sub("img.RData", "CoefBase.csv", arqs)

VAZMAX <- lapply(arqs, function(arq) {

    VAZMAX <- read.csv(arq)
    VAZMAX <- tail(VAZMAX[, 2], 1)
    cod    <-  read.table(sub("Niveis .*", "CodUsi.txt", arq))[1, 1]
    data.frame(USI = cod, VAZMAX = VAZMAX)
})

VAZMAX <- do.call(rbind, VAZMAX)
VAZMAX <- VAZMAX[order(VAZMAX[, 1]), ]

usethis::use_data(VAZMAX, internal = TRUE, overwrite = TRUE)
