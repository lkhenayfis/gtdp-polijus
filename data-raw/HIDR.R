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

usethis::use_data(HIDR, internal = TRUE, overwrite = TRUE)
