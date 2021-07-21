devtools::load_all()

##################################### DADOS INTERNOS DO PACOTE #####################################

# Preparo do hidr **COM CADASTRO ANTIGO** ----------------------------------------------------------

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

# VAZOES MAXIMAS PARA EXTENSAO ---------------------------------------------------------------------

arq <- "X:/Ciclo 2 - 2010 a 2019/Dados recebidos/_EPE/Vazoes Maximas Historicas 1931-2018.xlsx"
VAZMAX <- readxl::read_xlsx(arq)[, -2]
setDT(VAZMAX)
colnames(VAZMAX) <- c("USI", "VAZMAX")
setorder(VAZMAX, USI)


# ESCRITA DOS DADOS INTERNOS -----------------------------------------------------------------------

usethis::use_data(HIDR, VAZMAX, internal = TRUE, overwrite = TRUE)

##################################### DADOS EXTERNOS DO PACOTE #####################################

# EXEMPLO DE PLANILHA IMPORTADA --------------------------------------------------------------------

path <- system.file("extdata", "dummydata.xlsx", package = "polijus")

dummydata <- importadados(path)

usethis::use_data(dummydata, overwrite = TRUE)
