fasta_files = list.files(path="~tine/PRDM9_teleost/Teleostei_assemblies", pattern = "\\.fna$", recursive=TRUE, full.names=TRUE)
info_files = list.files(path="~tine/PRDM9_teleost/Teleostei_assemblies", pattern= "info.*assemblies", full.names=TRUE)

## because the 3rd one had to be edited manually
info_files[3] <- "info_teleostei_assemblies_ncbi_03_ed.tab"
ass_info_l <- lapply(info_files, read.delim, sep="\t", header=TRUE)

ass_info <- do.call( rbind, ass_info_l )

accessions <- sub(".+?/data/([^/]+).+", "\\1", fasta_files)
b <- nchar(accessions) != 15
accessions[b] <- sub(".+?_03/(GCA_[^_]+).+", "\\1", fasta_files[b])

## get rid of duplicated assemblies
fasta_files <- fasta_files[ !duplicated(accessions) ]
accessions <- accessions[ !duplicated(accessions) ]


##Remove cd files
names(fasta_files) <- accessions

fasta_files <- fasta_files[ !grepl("cds_from_genomic", fasta_files) ]
fasta_files <- fasta_files[ names(fasta_files) %in% ass_info$Assembly.Accession ]

dna_data <- cbind(ass_info[match(names(fasta_files), ass_info$Assembly.Accession), ], "file"=fasta_files)

write.csv(dna_data, "assemblies_info.csv")

