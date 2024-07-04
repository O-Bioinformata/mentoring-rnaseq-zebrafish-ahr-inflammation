# Instale pacotes se ainda não estiverem instalados
if (!require(digest)) install.packages("digest", dependencies = TRUE)
if (!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE)
if (!require(stringr)) install.packages("stringr", dependencies = TRUE)

# Carregar pacotes necessários
library(digest)
library(tidyverse)
library(stringr)

# Função para calcular o MD5 de um arquivo
calculate_md5 <- function(file_path) {
  md5 <- digest(file = file_path, algo = "md5")
  return(md5)
}

# Caminho do arquivo compactado e do CSV de saída
archive_path <- "D:/RNAseq_Raw data_Zebrafish_RRODRIGO_DISNER/part2/part2/"  # Substitua pelo caminho correto
output_csv <- "checksums_md5_part2.csv"

# Diretório temporário para extração
extract_to <- "temp_extracted"

# Cria o diretório de extração se não existir
if (!dir.exists(extract_to)) dir.create(extract_to)

# Lista todos os arquivos
file_list <- list.files(extract_to, pattern = "\\.fq\\.gz$", full.names = TRUE, recursive = TRUE)

# Calcula os MD5 e cria uma tibble com os resultados
md5_results <- tibble(
  Filename = basename(file_list),
  MD5_Checksum = sapply(file_list, calculate_md5)
)

# Salva os resultados em um arquivo CSV
write_csv(md5_results, output_csv)

