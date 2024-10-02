#Diana Escalante A01253468
#Libia Flores A01253425

#====Librerias=====
library(ape)
library(Biostrings)
library(DECIPHER)
library(seqinr)
library(ape)
library(ggtree)
library(viridis)
library(ggplot2)
library(ade4)

#====Directorio====
setwd("C:/Users/libia/OneDrive/Escritorio/2DOSEMESTRE/AnalisisdeBiologiaComputacional")
sars_cov2 <- c("MW973840","MW708309","MW731339","OM540761","OM540795","OM540759","OR856843","OR856840","OR856860","MZ312096","MW988204","MZ312090")

#====DNABin====
virus_seq <- read.GenBank(sars_cov2)
str(virus_seq)

#====Funcion write.dna====
write.dna(virus_seq, file = "SARS_COV2_seqs.fasta", format = "fasta")

#====Cargar secuencias en un nuevo vector====
sars_cov2_seq_not_align <- readDNAStringSet("SARS_COV2_seqs.fasta", format = "fasta")
sars_cov2_seq_not_align

#====Orientacion de nucleotidos====
sars_cov2_seq_not_align <- OrientNucleotides(sars_cov2_seq_not_align)

#====Alineamiento de las secuencias====
sars_cov2_seq_align <- AlignSeqs(sars_cov2_seq_not_align)

#====Visualizacion====
BrowseSeqs(sars_cov2_seq_align)

#====Guardar el resultado de alineamiento====
writeXStringSet(sars_cov2_seq_align, file = "SARS_COV2_seq_ALIGN.fasta")

#====Cargar archivo de alineamiento====
sars_cov_aligned <- read.alignment("SARS_COV2_seq_ALIGN.fasta", format = "fasta")

#====Crear matriz distancia====
matriz_distancia <- dist.alignment(sars_cov_aligned, matrix = "similarity")
matriz_distancia
temp <- as.data.frame(as.matrix(matriz_distancia))
temp

#====Tabla en escala de grises====
table.paint(temp, cleg=0, clabel.row=0.5, clabel.col=0.5) + scale_color_viridis()

#====Crear objeto phylo===
sars_cov_tree <- nj(matriz_distancia)

#====Crear arbol filogenetico====
sars_cov_tree <- ladderize(sars_cov_tree)
plot_sars_cov_filogenia <- ggtree(sars_cov_tree) + geom_tiplab() + ggtitle("Arbol filogenetico de los genomas SARS-COV2")
plot_sars_cov_filogenia


