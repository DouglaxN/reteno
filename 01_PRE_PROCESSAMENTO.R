######################################################################
#   Projeto: Análise de Expressão Diferencial (RNA-seq)
#   Contexto: Comparação entre grupo controle e reteno
#######################################################################

######################################################################
# 2️⃣ CRIAÇÃO DO DICIONÁRIO GENE ↔ TRANSCRITO
######################################################################

# Caminho do arquivo de anotação GTF comprimido
gtf_path <- "annotation/gencode.v47.annotation.gtf.gz"

# Caminho para o banco de dados SQLite que será criado a partir do GTF
txdb_sqlite <- "annotation/gencode.v47.annotation.gtf.sqlite"

# Se o banco SQLite ainda não existe, cria a partir do GTF
if (!file.exists(txdb_sqlite)) {
  # Cria o objeto TxDb a partir do GTF
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
  
  # Salva o banco de dados no formato SQLite
  saveDb(txdb, txdb_sqlite)
}

# Carrega o banco de dados SQLite para uso
txdb <- loadDb(txdb_sqlite)

# Extrai relação entre genes (GENEID) e transcritos (TXNAME)
tx_info <- AnnotationDbi::select(
  txdb,                   # Banco de dados carregado
  keys(txdb, "GENEID"),   # Lista de IDs de genes
  "TXNAME",               # Coluna desejada: nome dos transcritos
  "GENEID"                # Coluna-chave para mapeamento
)

# Conta o número de transcritos por gene e adiciona como coluna
tx_info$ntx <- table(tx_info$GENEID)[match(tx_info$GENEID, names(table(tx_info$GENEID)))]

# Cria dataframe no formato exigido pelo tximport: colunas "tx" e "gene"
tx2gene <- data.frame(
  tx = tx_info$TXNAME,    # Nome do transcrito
  gene = tx_info$GENEID,  # ID do gene associado
  stringsAsFactors = FALSE
)

# Salva o dicionário em formato RData para reutilização
save(tx2gene, file = "annotation/dicionario_gene_tx.RData")

######################################################################
# 3️⃣ METADADOS - GRUPO CONTROLE
######################################################################

# Lê metadados do grupo controle
controle_raw <- read_csv("metadado/meta_24h_controle.csv")

# Cria dataframe padronizado para o controle
metadado_controle <- data.frame(
  sample_id = paste0("Control_", 1:6),  # IDs únicos padronizados
  name      = controle_raw$sample,      # Nome original das amostras
  type      = "Control"                 # Define tipo = Controle
)

# Define sample_id como rownames (obrigatório para DESeq2)
rownames(metadado_controle) <- metadado_controle$sample_id

######################################################################
# 4️⃣ METADADOS - GRUPO OCCC
######################################################################

# Lê metadados do grupo OCCC
reteno_raw <- read_csv("metadado/meta_24h_reteno.csv")

# Cria dataframe padronizado para o controle
metadado_reteno <- data.frame(
  sample_id = paste0("Control_", 1:6),  # IDs únicos padronizados
  name      = controle_raw$sample,      # Nome original das amostras
  type      = "Reteno"                 # Define tipo = reteno
)

# Cria dataframe padronizado para o grupo OCCC
metadado_reteno <- data.frame(
  sample_id = reteno_raw$sample,               # IDs originais
  name      = substr(reteno_raw$fastq_1, 1, 20),# Extrai nome curto do FASTQ
  type      = "RETENO"                         # Define tipo = Câncer
)

# Define sample_id como rownames
rownames(metadado_occc) <- metadado_occc$sample_id

# Combina metadados controle + OCCC
metadado_final <- rbind(metadado_controle, metadado_occc)

# Salva metadados
save(metadado_final, file = "783540vs818977/RDatas/METADADO_OCCC.RData")

######################################################################
# 5️⃣ IMPORTAÇÃO E ORGANIZAÇÃO DAS TABELAS DE CONTAGEM
######################################################################

# Lê contagens de genes para OCCC e Controle
gene_counts_occc <- read.delim("dados_processados/dados_rnabulk_783540/analise_783540/tximport/gene_counts.tsv")
gene_counts_ctrl <- read.delim("dados_processados/dados_rnabulk_818977/tximport/gene_counts.tsv")

# Converte todos os valores para inteiros (necessário para DESeq2)
gene_counts_occc[] <- lapply(gene_counts_occc, as.integer)
gene_counts_ctrl[] <- lapply(gene_counts_ctrl, as.integer)

# Ajusta nomes das colunas para bater com os metadados
colnames(gene_counts_ctrl) <- metadado_controle$sample_id
colnames(gene_counts_occc) <- gsub("_T1$", "", colnames(gene_counts_occc))

# Ajusta colunas de OCCC para que coincidam com sample_id
colnames(gene_counts_occc)[match(metadado_occc$sample_id, colnames(gene_counts_occc))] <- metadado_occc$sample_id

# Adiciona coluna com ENSEMBL ID (necessário para merge)
gene_counts_occc <- rownames_to_column(gene_counts_occc, var = "ENSEMBLID")
gene_counts_ctrl <- rownames_to_column(gene_counts_ctrl, var = "ENSEMBLID")

# Junta contagens de Controle + OCCC
gene_counts_total <- merge(gene_counts_ctrl, gene_counts_occc, by = "ENSEMBLID")

# Remove genes duplicados
gene_counts_total <- gene_counts_total[!duplicated(gene_counts_total$ENSEMBLID), ]

# Define ENSEMBLID como rownames e remove coluna extra
rownames(gene_counts_total) <- gene_counts_total$ENSEMBLID
gene_counts_total$ENSEMBLID <- NULL

# Reordenando os metadados para bater com a ordem da matriz
metadado_final <- metadado_final[colnames(gene_counts_total), ]

# Salva tabela final de contagens
save(gene_counts_total, file = "783540vs818977/RDatas/GENE_COUNTS.RData")

######################################################################
# 6️⃣ ANÁLISE PCA COM DESeq2
######################################################################

# Converte a coluna type para fator (necessário para DESeq2)
metadado_final$type <- as.factor(metadado_final$type)

# Cria objeto DESeq2 com dados brutos
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(gene_counts_total), # Matriz de contagens
  colData   = metadado_final,               # Metadados
  design    = ~ type                        # Variável de interesse
)

# Normaliza dados usando VST (Variance Stabilizing Transformation)
dds_vst <- vst(dds)

# Extrai matriz de contagens normalizadas
counts_genes_vst <- assay(dds_vst)

# Garante ordem entre metadado e matriz
rownames(metadado_final) <- colnames(counts_genes_vst)
counts_genes_vst <- counts_genes_vst[, rownames(metadado_final)]

# Salva objetos para análises futuras
save(metadado_final, counts_genes_vst, dds, dds_vst,
     file = "783540vs818977/RDatas/DDS.RData")

######################################################################
# 7️⃣ VISUALIZAÇÕES (PCA, Screeplot, Dendrograma)
######################################################################

# Realiza PCA
pca_res <- pca(counts_genes_vst, metadata = metadado_final)

# Salva gráfico PCA
png("783540vs818977/imgs/PCA.png", width = 2500, height = 2000, res = 250)
biplot(pca_res, colby = "type", legendPosition = "bottom", title = "PCA - OCCC vs Controle")
dev.off()

# Salva gráfico Screeplot
png("783540vs818977/imgs/screeplot.png", width = 2500, height = 2000, res = 250)
screeplot(pca_res, title = "Screeplot - OCCC vs Controle")
dev.off()

# Prepara dados para clustering hierárquico
counts_dendo <- t(counts_genes_vst) # Transpõe para amostras como linhas
dist_matrix <- dist(counts_dendo)   # Calcula matriz de distâncias
cluster <- hclust(dist_matrix, method = "complete") # Agrupa por similaridade

# Salva dendrograma
png("783540vs818977/imgs/cluster.png", width = 3000, height = 2500, res = 300)
plot(cluster, main = "Cluster Hierárquico - Amostras")
dev.off()


