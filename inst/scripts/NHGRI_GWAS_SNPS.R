############################################################
###      Effects of NHGRI GWAS SNPs on branchpoints      ###
############################################################

###### download files ######

# download GTF and .fa files
system("wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz")
system("gunzip gencode.v21.annotation.gtf.gz")

# download NHGRI GWAS file
# date: 26/06/17
system("wget https://www.ebi.ac.uk/gwas/api/search/downloads/full")
system("mv full gwas_catalog_v1.0-associations_e89_r2017-06-19.tsv")
###### run branchpointer ######

library(biomaRt)
library(branchpointer)

options(stringsAsFactors = FALSE)
options(scipen = 999)

mart.snp <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp",host="www.ensembl.org")

# read in GWAS annotation
NHGRI_GWAS <- read.delim("gwas_catalog_v1.0-associations_e89_r2017-06-19.tsv",
                         stringsAsFactors = FALSE)

# manually create table -- biomaRt queries are slow & unreliable for large numbers
NHGRI_GWAS$allele <- NA
NHGRI_GWAS_filtered <- NHGRI_GWAS[-grep(";", NHGRI_GWAS$SNPS),]
NHGRI_GWAS_filtered <- NHGRI_GWAS_filtered[!duplicated(NHGRI_GWAS_filtered$SNPS),]


# queries of size 1000
split_size <- 1000
splits <- ceiling(nrow(NHGRI_GWAS_filtered) / split_size)

for(s in 1:splits){

  end_index <- s * split_size
  start_index <- end_index - (split_size - 1)
  end_index <- min(end_index, nrow(NHGRI_GWAS_filtered))
  
  snpInfo <- biomaRt::getBM(attributes = c("refsnp_id",'refsnp_source', "chr_name",
                                           "chrom_start", "allele"),
                            filters = "snp_filter", values = NHGRI_GWAS_filtered$SNPS[start_index:end_index], mart = mart.snp)

  m <- match(snpInfo$refsnp_id, NHGRI_GWAS_filtered$SNPS)
  NHGRI_GWAS_filtered$allele[m] <- snpInfo$allele
  
  message(s)
}

# check alleles -- only SNPs wanted
NHGRI_GWAS_filtered <- NHGRI_GWAS_filtered[which(!is.na(NHGRI_GWAS_filtered$allele)),]
NHGRI_GWAS_filtered$allele_ref <- unlist(lapply(stringr::str_split(NHGRI_GWAS_filtered$allele, "/"), "[[", 1))
NHGRI_GWAS_filtered <- NHGRI_GWAS_filtered[NHGRI_GWAS_filtered$allele_ref %in% c("A","C","G","T"),]
NHGRI_GWAS_filtered$allele_alt <- unlist(lapply(stringr::str_split(NHGRI_GWAS_filtered$allele, "/"), "[[", 2))

n_alt <- stringr::str_count(NHGRI_GWAS_filtered$allele, "/")
NHGRI_GWAS_filtered_copy <- NHGRI_GWAS_filtered[which(n_alt > 1),]
NHGRI_GWAS_filtered_copy$allele_alt <- unlist(lapply(stringr::str_split(NHGRI_GWAS_filtered_copy$allele, "/"), "[[", 3))
NHGRI_GWAS_filtered <- rbind(NHGRI_GWAS_filtered, NHGRI_GWAS_filtered_copy)
NHGRI_GWAS_filtered_copy <- NHGRI_GWAS_filtered[which(n_alt > 2),]
NHGRI_GWAS_filtered_copy$allele_alt <- unlist(lapply(stringr::str_split(NHGRI_GWAS_filtered_copy$allele, "/"), "[[", 4))
NHGRI_GWAS_filtered <- rbind(NHGRI_GWAS_filtered, NHGRI_GWAS_filtered_copy)
NHGRI_GWAS_filtered <- NHGRI_GWAS_filtered[NHGRI_GWAS_filtered$allele_alt %in% c("A","C","G","T"),]

NHGRI_GWAS_filtered$id <- with(NHGRI_GWAS_filtered, paste0(SNPS, ":", allele_ref, "/", allele_alt))
NHGRI_GWAS_filtered <- NHGRI_GWAS_filtered[which(NHGRI_GWAS_filtered$CHR_POS != ""),]

# write table for readQueryFile()
NHGRI_query_formatted <- data.frame(id = NHGRI_GWAS_filtered$id,
                                    chromosome = paste0("chr", NHGRI_GWAS_filtered$CHR_ID),
                                    chrom_start = NHGRI_GWAS_filtered$CHR_POS,
                                    strand = "x",
                                    ref_allele = NHGRI_GWAS_filtered$allele_ref,
                                    alt_allele = NHGRI_GWAS_filtered$allele_alt)

write.table(NHGRI_query_formatted, file="GWAS_SNP_query.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)

# format exon annotation
t1 <- system.time({exons <- gtfToExons("gencode.v21.annotation.gtf")})

t2 <- system.time({query <- readQueryFile("GWAS_SNP_query.txt", "SNP", exons)})
t2.noFilter <- system.time({query <- readQueryFile("GWAS_SNP_query.txt", "SNP", exons, filter = F)})

# predict branchpoints
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
t3 <- system.time({branchpointPredictions <- predictBranchpoints(query,queryType = "SNP",BSgenome = g)})

# filter for SNPs that create or delete branchpoints
t4 <- system.time({summary <- predictionsToSummary(query,branchpointPredictions)})

