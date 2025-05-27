#!/usr/bin/env Rscript

# MTP-DB is our transportome-wide SQLite database, assembled by Daedalus and
# used by Ariadne to build up all the gene lists of interest.

# --- Packages -----------------------------------------------------------------

library(dplyr)
library(DBI)
library(RSQLite)

# --- Functions ----------------------------------------------------------------

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments
mtpdb <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# live test (from the project root directory)
mtpdb <- "./data/MTPDB.sqlite"
out_dir <- "./data/out"

if (Sys.info()["sysname"] == "Windows") {
  # Can't access the DB directly on WSL when running from Windows...
  mtpdb <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "MTPDB.sqlite")
  file.copy(from = "./data/MTPDB.sqlite", to = mtpdb)
}





connection <- dbConnect(SQLite(), dbname = mtpdb)



query_pores <- 
  "SELECT DISTINCT
  	channels.ensg,
  	gene_names.hugo_gene_symbol,
  	gene_names.hugo_gene_name
  FROM
  	channels
  	JOIN
  		gene_names
  		ON channels.ensg = gene_names.ensg
  UNION SELECT DISTINCT
  	aquaporins.ensg,
  	gene_names.hugo_gene_symbol,
  	gene_names.hugo_gene_name
  FROM
  	aquaporins
  	JOIN
  		gene_names
  		ON aquaporins.ensg = gene_names.ensg
  ORDER BY gene_names.hugo_gene_symbol;"


query_pumps <-
  "SELECT DISTINCT
  	pumps.ensg,
  	gene_names.hugo_gene_symbol,
  	gene_names.hugo_gene_name
  FROM
  	pumps
  	JOIN
  		gene_names
  		ON pumps.ensg = gene_names.ensg;"


query_iSLCs <- 
  "SELECT DISTINCT
  	solute_carriers.ensg,
  	solute_carriers.carried_solute,
  	gene_names.hugo_gene_symbol,
  	gene_names.hugo_gene_name
  FROM
  	solute_carriers
  	JOIN gene_names
  		ON solute_carriers.ensg = gene_names.ensg
  WHERE NOT EXISTS (
  	SELECT 1
  	FROM solute_carriers AS sc_sub
  	WHERE sc_sub.ensg = solute_carriers.ensg
  	  AND sc_sub.carried_solute NOT IN ('ion', 'cation', 'anion', 'Na+', 'K+', 'Ca2+', 'Cl-', 'H+', 'HCO3-', 'OH-', 'Fe2+', 'Mg2+', 'Zn2+', 'Mn2+', 'Li+')
  );"




pores <- dbGetQuery(connection, query_pores)

iSLCs <- dbGetQuery(connection, query_iSLCs)

write.csv(pores, row.names = FALSE, file = file.path(out_dir, "ICT_pores.csv"))














dbDisconnect(connection)





