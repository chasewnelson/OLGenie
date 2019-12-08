############################################################################################################
### Supplementary Material for Nelson CW, Ardern Z, Wei X, "OLGenie: Estimating Natural Selection to Predict Functional Overlapping Genes"
### Convert XML alignment from BLAST to FASTA

library(tidyverse)
library(GenomicFeatures)
library(Biostrings)
library(xml2)

# LOAD
input_file_path <- "ENTER/YOUR/PATH/HERE"
my_XML_file <- xml2::read_xml(x = input_file_path)

Hit_id <- xml_find_all(my_XML_file, "//Hit_id")
Hit_id_char <- as.character(Hit_id)
Hit_id_char <- str_replace_all(string = Hit_id_char, pattern = "<Hit_id>", replacement = "")
Hit_id_char <- str_replace_all(string = Hit_id_char, pattern = "</Hit_id>", replacement = "")

Hit_def <- xml_find_all(my_XML_file, "//Hit_def")
Hit_def_char <- as.character(Hit_def)
Hit_def_char <- str_replace_all(string = Hit_def_char, pattern = "<Hit_def>", replacement = "")
Hit_def_char <- str_replace_all(string = Hit_def_char, pattern = "</Hit_def>", replacement = "")
Hit_def_char_short <- str_sub(string = Hit_def_char, start = 1, end = 100)

Hit_accession <- xml_find_all(my_XML_file, "//Hit_accession")
Hit_accession_char <- as.character(Hit_accession)
Hit_accession_char <- str_replace_all(string = Hit_accession_char, pattern = "<Hit_accession>", replacement = "")
Hit_accession_char <- str_replace_all(string = Hit_accession_char, pattern = "</Hit_accession>", replacement = "")

Hit_len <- xml_find_all(my_XML_file, "//Hit_len")
Hit_len_char <- as.character(Hit_len)
Hit_len_char <- str_replace_all(string = Hit_len_char, pattern = "<Hit_len>", replacement = "")
Hit_len_char <- str_replace_all(string = Hit_len_char, pattern = "</Hit_len>", replacement = "")
Hit_len_int <- as.integer(Hit_len_char)

Hsp_num <- xml_find_all(my_XML_file, "//Hsp_num")
Hsp_num_char <- as.character(Hsp_num)
Hsp_num_char <- str_replace_all(string = Hsp_num_char, pattern = "<Hsp_num>", replacement = "")
Hsp_num_char <- str_replace_all(string = Hsp_num_char, pattern = "</Hsp_num>", replacement = "")
Hsp_num_int <- as.integer(Hsp_num_char)

Hsp_align_len <- xml_find_all(my_XML_file, "//Hsp_align-len")
Hsp_align_len_char <- as.character(Hsp_align_len)
Hsp_align_len_char <- str_replace_all(string = Hsp_align_len_char, pattern = "<Hsp_align-len>", replacement = "")
Hsp_align_len_char <- str_replace_all(string = Hsp_align_len_char, pattern = "</Hsp_align-len>", replacement = "")
Hsp_align_len_int <- as.integer(Hsp_align_len_char)
Hsp_align_len_int_Hit1 <- Hsp_align_len_int[Hsp_num_int == 1]

Hsp_hseq <- xml_find_all(my_XML_file, "//Hsp_hseq")
Hsp_hseq_char <- as.character(Hsp_hseq)
Hsp_hseq_char <- str_replace_all(string = Hsp_hseq_char, pattern = "<Hsp_hseq>", replacement = "")
Hsp_hseq_char <- str_replace_all(string = Hsp_hseq_char, pattern = "</Hsp_hseq>", replacement = "")

HEADERS <- paste0(Hit_id_char, ' ', Hit_accession_char, ' ', Hit_def_char_short, ' length=', Hit_len_int, ' aln_length=', Hsp_align_len_int_Hit1)
HEADERS <- str_replace_all(string = HEADERS, pattern = "\\|", replacement = "_")
HEADERS <- str_replace_all(string = HEADERS, pattern = "_ ", replacement = " ")

# Construct DNA String Set: only select the first hits
Hsp_hseq_char_Hit1 <- Hsp_hseq_char[Hsp_num_int == 1]
Hits_DSS <- DNAStringSet(x = Hsp_hseq_char_Hit1)
names(Hits_DSS) <- HEADERS

# Define output fasta
output_file_path <- input_file_path
output_file_path <- str_replace(string = output_file_path, pattern = ".xml", replacement = ".fasta")
#writeXStringSet(x = Hits_DSS, filepath = output_file_path, format = 'fasta')

# NOW USE process_alignment.R

