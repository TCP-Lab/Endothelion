


# Validation workflow

verficare la presenza dei prefissi dentro la cartella Counts:	tree Counts/
verificare la presenza del prefisso nel file di CountMatrix_genes_TPM.tsv
verficare la presenza dei prefissi dentro la cartella MultiQC_out:	tree MultiQC_out/

aggiungere PCA report: qcfastq --tool=PCA .
per il controllo di qualità scaricare:
		il MultiQC report
		PCA report
scp fear@130.192.101.241:/home/fear/WORKS/Endothelion/Lines/hCMEC_D3/GSE138309/MultiQC_out/GSE138309_multiqc_report.html .

scp fear@130.192.101.241:/home/fear/WORKS/Endothelion/Lines/hCMEC_D3/GSE138309/PCA_out/GSE138309_CountMatrix_genes_TPM/* .

riscaricare i metadati e editare il campo extra da nano (0 per le run da escludere)
	 metaharvest -x=1 -e GSE138309 > GSE138309_meta.csv

cancella eventuali FASTQ file ancora presenti
	rm *.fastq.gz

crea l'archivio tar.gz dalla parent directory
	tar -czvf ./GSE205739/GSE205739.tar.gz GSE205739

scarica in locale i file da mandare a Zenodo
	scp fear@130.192.101.241:/home/fear/WORKS/Endothelion/Lines/hCMEC_D3/GSE205739/GSE205739* .

# Docs



SRA @ https://www.ncbi.nlm.nih.gov/sra
search: hCMEC D3
BioProject    12


OK - in now
PRJNA307652			GSE76528
PRJNA575504			GSE138309
PRJNA578611			GSE139133
PRJNA777606			GSE187565
PRJNA802135			GSE195781
PRJNA847413			GSE205739
PRJEB48614			E-MTAB-11129
PRJNA667281			--
PRJNA896725			--

excluded
PRJNA607654			GSE145581		miRNA-Seq
PRJNA1073892		GSE255171		regulatory T cells (Tregs)
PRJNA307651			GSE76530		miRNA-Seq





Nello studio PRJNA777606 (GSE187565) 2 run (SRR16764424 e SRR16764425) sono in effetti
bulk RNA-Seq, mentre 2 (SRR16764426 e SRR16764427) erano invece ATAC-Seq erroneamente incluse (ora rimosse).


---------------
GSE205739
allinea intorno al 50%, ma tutto il resto sembra ok
---------------

Of the 5 control runs of PRJNA802135 / GSE195781 (2022) 3 turned out to be
identical to runs already published in the previous study of the same group
of authors (2016) with ID PRJNA307652 / GSE76528

PRJNA802135			PRJNA307652
(GSE195781)			(GSE76528)
2022				2016

SRR17833475		=	SRR3085451
SRR17833476		=	SRR3085449
SRR17833478		=	SRR3085446

risultano effettivamente identiche le reads all'interno dei rispettivi FASTQ,
e le diverse dimensioni dipendono solo dal diverso pattern usato per la linea di
intestazione di ogni read
zcat SRR17833478_1.fastq.gz | wc -l
zcat SRR3085446_1.fastq.gz | wc -l

zcat SRR17833478_1.fastq.gz | head
zcat SRR3085446_1.fastq.gz | head

zcat SRR17833478_1.fastq.gz | tail
zcat SRR3085446_1.fastq.gz | tail

questo spiega l'effetto di batch che separa quelle 3 run dagli altri 2 controlli
gli anomali livelli di deviazione standard che si producono mettendo tutti i
sample insieme e l'improbabile livello di correlazione tra i 2 dataset nel
momento in cui si eliminano le 2 run diverse.

> model_mean$Std_Dev |> mean()
[1] 0.1720537

model |> lapply(geneStats) -> series_mean

series_mean |> sapply(\(x) mean(x$Std_Dev))
 GSE138309  GSE139133  GSE195781  GSE205739   GSE76528 
0.09534965 0.05233470 0.06666573 0.09942313 0.11305201 


> model |> geneStats() |> filter(Mean > 1) -> model_mean
> mean(model_mean$Std_Dev)
0.5259145

> model |> lapply(geneStats) |> lapply(filter, Mean > 1) -> series_mean
> series_mean |> sapply(\(x) mean(x$Std_Dev))
GSE138309 GSE139133 GSE195781 GSE205739  GSE76528 
0.2290709 0.1265768 0.1474828 0.2312040 0.3172325 


---------
GSE195781 è da escludere completamente perché
3 campioni su 5 (SRR17833475, SRR17833476, SRR17833478) appartengono in realtà ad un altro studio
mentre i 2 rimanenti hanno gravi problemi di contaminazione da adapter, con % di allineamento del 20%
---------


