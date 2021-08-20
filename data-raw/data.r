devtools::load_all()

##################################### DADOS INTERNOS DO PACOTE #####################################

# Preparo do hidr **COM CADASTRO ANTIGO** ----------------------------------------------------------

HIDR <- read.csv("data-raw/HIDR_PMO_Jun2021.csv", fill = TRUE, sep = ";", dec = ",", stringsAsFactors = FALSE)

# Identifica quais colunas contem numero de maquinas por grupo
colnummaq <- grep("X\\.Maq\\.[[:digit:]]\\.", colnames(HIDR))

# Identifica as colunas que contem a vazao efetiva de cada grupo
colvazef <- grep("QEf\\.[[:digit:]]\\.", colnames(HIDR))

# Vazao efetiva da usina
vazef <- rowSums(HIDR[, colnummaq] * HIDR[, colvazef])

# Namax da usina a jusante
namaxjus <- as.numeric(sub(" - .*", "", HIDR[, "Jusante"]))
namaxjus <- match(namaxjus, HIDR[, 1])
namaxjus <- HIDR[namaxjus, 11]
namaxjus[namaxjus == 0] <- NA_real_

# Identifica as colunas que contem os coeficientes dos polinomios de jusante
colpvn <- grep("PJA[[:digit:]]\\.[[:digit:]]", colnames(HIDR))

HIDR <- cbind(USI = HIDR[, 1], VAZEF = vazef, NAMAX = namaxjus, HIDR[, colpvn])

dummy <- HIDR[1, ]
dummy$USI <- 999
dummy$VAZEF <- 440
dummy$NAMAX <- 34
dummy[, 4] <- 34
dummy[, 5] <- 1 / 1400

HIDR <- rbind(HIDR, dummy)
setDT(HIDR)

# VAZOES MAXIMAS PARA EXTENSAO ---------------------------------------------------------------------

arq <- "X:/Ciclo 2 - 2010 a 2019/Dados recebidos/_EPE/Vazoes Maximas Historicas 1931-2018.xlsx"
VAZMAX <- readxl::read_xlsx(arq)[, -2]
setDT(VAZMAX)
colnames(VAZMAX) <- c("USI", "VAZMAX")
setorder(VAZMAX, USI)


# ESCRITA DOS DADOS INTERNOS -----------------------------------------------------------------------

usethis::use_data(HIDR, VAZMAX, internal = TRUE, overwrite = TRUE)

##################################### DADOS EXTERNOS DO PACOTE #####################################

# EXEMPLO DE PLANILHA IMPORTADA E TRATADA ----------------------------------------------------------

path <- system.file("extdata", "dummydata.xlsx", package = "polijus")

dummydata <- importadados(path)
dummydata <- filtravazest(dummydata)
dummydata <- classfiltrapats(dummydata)
dummydata <- evalremanso(dummydata)

usethis::use_data(dummydata, overwrite = TRUE)
