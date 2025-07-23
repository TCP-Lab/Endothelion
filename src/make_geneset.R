#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

library(dplyr, warn.conflicts = FALSE)
library(DBI)
library(RSQLite)

# --- Functions ----------------------------------------------------------------

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments
db_path <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# # Live debug (from the project root directory)
# db_path <- "./data/MTPDB.sqlite"
# out_dir <- "./data/out"

if (Sys.info()["sysname"] == "Windows") {
  # Can't access the DB directly on WSL when running from Windows...
  db_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "MTPDB.sqlite")
  file.copy(from = "./data/MTPDB.sqlite", to = db_path)
}

# Connect to the MTP-DB
connection <- dbConnect(SQLite(), dbname = db_path)

# Pores (Ion Channels + AQPs)
query_pores <-
    "SELECT DISTINCT
    	channels.ensg,
    	gene_names.hugo_gene_symbol,
    	gene_names.hugo_gene_name
    FROM
    	channels JOIN gene_names ON channels.ensg = gene_names.ensg
    
    UNION
    
    SELECT DISTINCT
    	aquaporins.ensg,
    	gene_names.hugo_gene_symbol,
    	gene_names.hugo_gene_name
    FROM
    	aquaporins JOIN gene_names ON aquaporins.ensg = gene_names.ensg
    
    ORDER BY gene_names.hugo_gene_symbol"

# Pumps
query_pumps <-
    "SELECT DISTINCT
      	pumps.ensg,
      	gene_names.hugo_gene_symbol,
      	gene_names.hugo_gene_name
    FROM
        pumps JOIN gene_names ON pumps.ensg = gene_names.ensg
    
    ORDER BY gene_names.hugo_gene_symbol"

# Inorganic SLCs
inorganics <- paste(
    "'ion", # mind the opening single quote !!
    "cation",
    "anion",
    "phosphate",
    "Na+",
    "K+",
    "Ca2+",
    "Cl-",
    "H+",
    "HCO3-",
    "OH-",
    "Fe2+",
    "Fe3+",
    "I-",
    "H2PO4-",
    "HPO42-",
    "SO42-",
    "NH4+",
    "Mg2+",
    "Zn2+",
    "Mn2+",
    "Li+",
    "cisplatin",
    "copper",
    "Cd2+'", # mind the terminal single quote !!
    sep = "', '")

query_iSLCs <- paste0(
    "SELECT DISTINCT
    	solute_carriers.ensg,
    	solute_carriers.carried_solute,
      	gene_names.hugo_gene_symbol,
      	gene_names.hugo_gene_name
    FROM
    	solute_carriers JOIN gene_names ON solute_carriers.ensg = gene_names.ensg
    WHERE
    	solute_carriers.ensg NOT IN (
    		SELECT DISTINCT
    			ensg
    		FROM
    			solute_carriers
    		WHERE
    			carried_solute NOT IN (", inorganics, ")
    		)
    	AND solute_carriers.carried_solute IS NOT NULL
    
    ORDER BY hugo_gene_symbol")

# # Mr.Hedmad's solution
# query_iSLCs_alternative <- paste0(
#     "SELECT DISTINCT
#     	solute_carriers.ensg,
#     	solute_carriers.carried_solute,
#       	gene_names.hugo_gene_symbol,
#       	gene_names.hugo_gene_name
#     FROM
#     	solute_carriers JOIN gene_names ON solute_carriers.ensg = gene_names.ensg
#     WHERE
#     	solute_carriers.ensg IN (
#     		SELECT DISTINCT
#     			ensg
#     		FROM
#     			solute_carriers
#     		WHERE
#     			carried_solute IN (", inorganics, ")
#     		
#     		EXCEPT
#     		
#     		SELECT DISTINCT
#     			ensg
#     		FROM
#     			solute_carriers
#     		WHERE
#     			carried_solute IN (
#     				SELECT DISTINCT
#     					carried_solute
#     				FROM
#     					solute_carriers
#     				WHERE
#     					carried_solute NOT IN (", inorganics, ")
#     		)
#     )")

# Make the calls
pores <- dbGetQuery(connection, query_pores)
pumps <- dbGetQuery(connection, query_pumps)
iSLCs <- dbGetQuery(connection, query_iSLCs)

# Temporary patch for iSLCs and ICs as retrieved from MTP-DB ver. 1.25.24
patch <- c(
    paste0("SLC30A", c(1:10)),
    paste0("SLC4A", c(1:11)),
    "SLC24A1",
    "CATSPERG",
    "CATSPERB",
    "CATSPERD",
    "CATSPERE",
    "CATSPERZ")

# Receptor Of Interest (ROIs)
ROIs <- c(
    "ADORA1",
    "ADORA2A",
    "ADORA2B",
    "ADORA3",
    "ADRA1A",
    "ADRA1B",
    "ADRA1D",
    "ADRA2A",
    "ADRA2B",
    "ADRA2C",
    "ADRB1",
    "ADRB2",
    "ADRB3",
    "CHRM1",
    "CHRM2",
    "CHRM3",
    "CHRM4",
    "CHRM5",
    "EGFR",
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "FGFR4",
    "FLT1",
    "FLT3",
    "FLT4",
    "GABBR1",
    "GABBR2",
    "GRM1",
    "GRM2",
    "GRM3",
    "GRM4",
    "GRM5",
    "GRM6",
    "GRM7",
    "GRM8",
    "HRH1",
    "HRH2",
    "HRH3",
    "HRH4",
    "IGF1R",
    "MAGT1",
    "P2RY1",
    "P2RY11",
    "P2RY12",
    "P2RY13",
    "P2RY14",
    "P2RY2",
    "P2RY4",
    "P2RY6",
    "PDGFRA",
    "PDGFRB")

# Extract Gene Symbols, combine, and save
c(pores$hugo_gene_symbol,
  pumps$hugo_gene_symbol,
  iSLCs$hugo_gene_symbol, patch, ROIs) |>
    unique() |> na.omit() |>
    write.table(sep = ",",
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE,
                file = file.path(out_dir, "ICT_set_v2.csv"))

# Disconnect from the MTP-DB
dbDisconnect(connection)

# 
# 
# c(pores$hugo_gene_symbol, patch) |>
#     unique() |> na.omit() |>
#     write.table(sep = ",",
#                 col.names = FALSE,
#                 row.names = FALSE,
#                 quote = FALSE,
#                 file = "tempo_ICs.csv")
# 


