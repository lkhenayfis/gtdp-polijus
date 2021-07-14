devtools::load_all()

########################## Preparo do hidr **COM CADASTRO ANTIGO**

HIDR <- read.csv("data-raw/HIDR.csv", fill = TRUE, sep = ";", dec = ",", stringsAsFactors = FALSE)

# Identifica quais colunas contem numero de maquinas por grupo
colnummaq <- grep("X\\.Maq\\.[[:digit:]]\\.", colnames(HIDR))

# Identifica as colunas que contem a vazao efetiva de cada grupo
colvazef <- grep("QEf\\.[[:digit:]]\\.", colnames(HIDR))

# Vazao efetiva da usina
vazef <- rowSums(HIDR[, colnummaq] * HIDR[, colvazef])

# Identifica as colunas que contem os coeficientes dos polinomios de jusante
colpvn <- grep("PJA[[:digit:]]\\.[[:digit:]]", colnames(HIDR))

HIDR <- cbind(USI = HIDR[, 1], VAZEF = vazef, HIDR[, colpvn])

dummy <- HIDR[1, ]
dummy$USI <- 999
dummy$VAZEF <- 950
dummy[, 3] <- 9
dummy[, 4] <- 1 / 1400

HIDR <- rbind(HIDR, dummy)

########################## VAZOES MAXIMAS PARA AJUSTE 

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

usethis::use_data(HIDR, VAZMAX, internal = TRUE, overwrite = TRUE)

########################## EXEMPLO DE DADO IMPORTADO

path <- system.file("extdata", "dummydata.xlsx", package = "polijus")

dummydata <- importadados(path)

usethis::use_data(dummydata, overwrite = TRUE)
