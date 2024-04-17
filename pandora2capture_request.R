#!/usr/bin/env Rscript

## Need sidora.core to pull library metainfo
if (!require('sidora.core')) {
  if(!require('remotes')) install.packages('remotes')
  remotes::install_github('sidora-tools/sidora.core', quiet=T)
} else {library(sidora.core)}

require(purrr)
require(tidyverse, warn.conflicts = F)
require(optparse)
require(readr)
require(stringr)

validate_csv <- function(csv) {
  probe_set <- c("1240k","16S-rRNA","Archaic Admixture","Big Yoruba","Borrelia/Borreliella","Brucella","Candida albicans",
                 "Canid","Capra (goat)","Cowpox","Fermentation (bacteria)","Fermentation (yeast)","Helicobacter pylori",
                 "Hepatitis B","Human Screening","Immunocapture","Mitochondrial","Mitochondrial Bison",
                 "Mitochondrial_cow (Tübingen)","Mitochondrial_deer (Tübingen)","Mix Mitochondrial/Human Screening",
                 "Mycobacterium leprae","Mycobacterium tuberculosis","oralStrep-ABP","oralStrep-genomes","Ovis (sheep)",
                 "Parvovirus B19","Plasmodium Mitogenome","Plasmodium Nuclear","Plasmodium Plastid (Apicoplast)","Salmonella enterica",
                 "Staphylococcus aureus","Streptococcus mutans","Streptococcus pyogenes","Treponema pallidum",
                 "Treponema succinifaciens","Y-Chromosom","Yersinia pestis")
  labs <- c("CoreUnit","JenaLab")
  sequencing_kits <- c("HiSeq4000 75 SR","HiSeq4000 75 PE**","HiSeq4000 150 PE**",
                     "Flex 75 PE*", "Flex 75 SR*",
                     "MiSeqNano 150 PE", 
                     "MiSeqMicro 150 PE","MiSeq 75 PE", "MiSeq 150 PE","MiSeq 300 PE", 
                     "NextSeq500MO 75 PE", "NextSeq500MO 75 SR", 
                     "NextSeq500MO 150 SR", "NextSeq500MO 150 PE**", 
                     "NextSeq500HO 75 PE**", "NextSeq500HO 75 SR**",
                     "NextSeq500HO 150 PE**",
                     "NovaSeqX 75 PE*",
                     "NovaSeqX 100PE*")
  sequencing_depths <- c("1M","2M","5M","10M","20M", "40M")
  incorrect_probe <- csv %>% filter(!Probe_Set %in% probe_set) %>% pull(Probe_Set)
  incorrect_lab <- csv %>% filter(!Lab %in% labs) %>% pull(Lab)
  incorrect_kit <- csv %>% filter(!Sequencing_Kit %in% sequencing_kits) %>% pull(Sequencing_Kit)
  incorrect_depth <- csv %>% filter(!Sequencing_Depth %in% sequencing_depths) %>% pull(Sequencing_Depth)
  
  if (length(incorrect_probe) == 0 & length(incorrect_lab) == 0 & length(incorrect_kit) == 0 & length(incorrect_depth) == 0) {
    cat("\n[pandora2capture_request.R]: CSV check passed")
  } else {
    if ( length(incorrect_probe) != 0 ) {
      cat(paste("\n[pandora2capture_request.R] error: Invalid probe set: '", incorrect_probe, 
           "'\nAccepted values: ", paste(probe_set,collapse=", "),"\n\n"))
    } 
    
   if ( length(incorrect_lab) != 0 ) {
      cat(paste("\n[pandora2capture_request.R] error: Invalid lab: '", incorrect_lab, 
           "'\nAccepted values: ", paste(labs,collapse=", "), "\n\n", sep = ""))
   }
   if ( length(incorrect_kit) != 0 ) {
      cat(paste("\n[pandora2capture_request.R] ERROR - Invalid sequencing kit: '", incorrect_kit, 
           "'\nAccepted values: ", paste(sequencing_kits,collapse=", "),"\n\n",sep = ""))
   }
   if ( length(incorrect_depth) != 0 ) {
      cat(paste("\n[pandora2capture_request.R] error: Invalid sequencing depth: '", incorrect_depth, 
           "'\nAccepted values: ", paste(sequencing_depths,collapse=", "),"\n\n"))
   }
   if ( length(incorrect_probe) != 0 | length(incorrect_lab) != 0 | length(incorrect_kit) != 0 | length(incorrect_depth) != 0 ) {
     stop()
   }
  }
}

validate_probe_set <- function(csv) {
  valid_entries <- c("1240k","16S-rRNA","Archaic Admixture","Big Yoruba","Borrelia/Borreliella","Brucella","Candida albicans",
                     "Canid","Capra (goat)","Cowpox","Fermentation (bacteria)","Fermentation (yeast)","Helicobacter pylori",
                     "Hepatitis B","Human Screening","Immunocapture","Mitochondrial","Mitochondrial Bison",
                     "Mitochondrial_cow (Tübingen)","Mitochondrial_deer (Tübingen)","Mix Mitochondrial/Human Screening",
                     "Mycobacterium leprae","Mycobacterium tuberculosis","oralStrep-ABP","oralStrep-genomes","Ovis (sheep)",
                     "Parvovirus B19","Plasmodium Mitogenome","Plasmodium Nuclear","Plasmodium Plastid (Apicoplast)","Salmonella enterica",
                     "Staphylococcus aureus","Streptococcus mutans","Streptococcus pyogenes","Treponema pallidum",
                     "Treponema succinifaciens","Y-Chromosom","Yersinia pestis")
  incorrect <- csv %>% filter(!probe_set %in% valid_entries)
  ifelse(nrow(incorrect) > 0, stop(call.=F, "\n[pandora2capture_request.R] error: Invalid sequencing kit: '", incorrect$Probe_Set, 
                                   "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"), incorrect)
  }

validate_lab <- function(value) {
  valid_entries <- c("CoreUnit","JenaLab") 
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[pandora2capture_request.R] error: Invalid lab: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

validate_kit <- function(csv) {
  valid_entries <- c("HiSeq4000 75 SR","HiSeq4000 75 PE**","HiSeq4000 150 PE**",
                     "Flex 75 PE*", 
                     "Flex 75 SR*",
                     "MiSeqNano 150 PE", 
                     "MiSeqMicro 150 PE","MiSeq 75 PE", "MiSeq 150 PE","MiSeq 300 PE", 
                     "NextSeq500MO 75 PE", "NextSeq500MO 75 SR", 
                     "NextSeq500MO 150 SR", "NextSeq500MO 150 PE**", 
                     "NextSeq500HO 75 PE**", "NextSeq500HO 75 SR**",
                     "NextSeq500HO 150 PE**",
                     "NovaSeqX 75 PE*",
                     "NovaSeqX 100PE*")
  incorrect <- csv %>% filter(!Sequencing_Kit %in% valid_entries)
  ifelse(nrow(incorrect) > 0, stop(call.=F, "\n[pandora2capture_request.R] error: Invalid sequencing kit: '", incorrect$Sequencing_Kit, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"), return(incorrect))
}

validate_depth <- function(value) { 
  valid_entries <- c("1M","2M","5M","10M","20M", "40M")
  incorrect <- csv %>% filter(!Sequencing_Depth %in% valid_entries)
  ifelse(nrow(incorrect) > 0, stop(call.=F, "\n[pandora2capture_request.R] error: Invalid sequencing kit: '", incorrect$Sequencing_Kit, 
                                   "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"), return(incorrect))
  
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[pandora2capture_request.R] error: Invalid depth: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

validate_libraries <- function(option, opt_str, value, parser) {
  ifelse(!is.na(value), 
         return(value), 
         stop(call.=F, "\n[pandora2capture_request.R] error: No file with libraries provided!\n Please provide a file with library IDs using the -f option"))
}

## Parse arguments ----------------------------
parser <- OptionParser(usage = "%prog [options] .credentials")
##Make -f mandatory with valid_entries
parser <- add_option(parser, c("-f", "--file_libraries_id"), type = 'character', 
                     action = "callback", dest = "file_libraries_id",
                     callback = validate_libraries,
                     help = "CSV file containing all the libraries ID, probe set, sequecing kit, sequencing depth, cost center and lab",
                     default = NA)
parser <- add_option(parser, c("-r", "--researcher"), type = 'character',
                     action = "store", dest = "researcher", default=NA,
                     help = "Name of the contact person for the request"
                     )
parser <- add_option(parser, c("-e", "--emailresearcher"), type = 'character',
                     action = "store", dest = "email", default=NA,
                     help = "Email of the contact person for the request"
                     )
parser <- add_option(parser, c("-o", "--outDir"), type = 'character',
                     action = "store", dest = "outdir",
                     help= "The desired output directory. By default, it is the current directory.",
                     default = "."
                     )
parser <- add_option(parser, c("-v", "--version"), action= 'store_true', default=TRUE, help = "Prints version")

arguments <- parse_args(parser, positional_arguments = 1)
opts <- arguments$options

if( opts$version == TRUE ) {
  print("[pandora2capture_request.R] version 1.0")
  quit()
}

cred_file <- arguments$args
libraries <- read.csv(opts$file_libraries_id, header = T)

#libraries

#Check all columns contain correct information: 
validate_csv(libraries)

#probe_set <- opts$probe_set
#lab <- opts$lab
contact <- paste(opts$researcher, opts$email, sep = ",")
#kit <- opts$kit
#depth <- opts$depth
#costCenter <- opts$costCenter
requestDate <- Sys.Date()

con <- sidora.core::get_pandora_connection(cred_file)



libraryTab <- get_df("TAB_Library", con) %>%
  convert_all_ids_to_values(con) %>%
  filter(library.Full_Library_Id %in% libraries$Library_ID)

JenaLab <- libraries %>%
  filter(Lab == "JenaLab")
  
  
CoreUnit <- libraries %>%
  filter(Lab == "CoreUnit")

if (nrow(JenaLab) > 0) {
  JenaTSV <- libraryTab %>%
    filter(library.Full_Library_Id %in% JenaLab$Library_ID) %>%
    left_join(JenaLab, by = c("library.Full_Library_Id" = "Library_ID")) %>%
    mutate(
      Library_ID = library.Full_Library_Id, 
      Researchers_Name = contact,
      Project_Name = ifelse(is.na(Cost_Center), library.Projects, as.character(Cost_Center)),
      Batch = library.Batch,
      Robot = case_when(
        library.Robot == FALSE ~ "n",
        library.Robot == TRUE ~ "y",
        is.na(library.Robot) ~ "NA"
      ),
      Protocol = library.Protocol,
      Probe_Set = Probe_Set,
      Position_on_plate = library.Position_on_Plate,
      Location = library.Location_Room,
      Drawer = paste(" ",library.Location, sep = ""),
      P7_Index_Name = library.P7_Index_Id,
      P7_Index_Sequence = library.P7_Index_Sequence,
      P5_Index_Name = library.P5_Index_Id,
      P5_Index_Sequence = library.P5_Index_Sequence,
      P5_Barcode = library.P5_Barcode_Sequence,
      P7_Barcode =library.P7_Barcode_Sequence,
      concentration_from_Nanodrop = "",
      additional_notes = "",
           ) %>%
    select(!contains("library.")) %>%
    select(!contains(c("Sequencing_Kit","Sequencing_Depth","Cost_Center","Lab", "MPI_Probe_Pool"))) %>%
    relocate(Probe_Set, .after = Protocol)
  nameOutput <- paste(paste("CaptureRequest","JenaLab",requestDate, sep = "_"),"tsv",sep = ".")
  outputTSV <- paste(opts$outdir,nameOutput, sep = "/")
  write_tsv(JenaTSV, file = outputTSV)
} 
if (nrow(CoreUnit) > 0) {
  Probe_Sets_Core_Unit <- read.csv("Probe_Set.csv") %>%
    select(Name,MPI_Probe_Pool)
  code_probe_set <- Probe_Sets_Core_Unit %>%
    filter(Name %in% CoreUnit$Probe_Set)
  CoreUnitTSV <- libraryTab %>%
    filter(library.Full_Library_Id %in% CoreUnit$Library_ID) %>%
    left_join(CoreUnit, by = c("library.Full_Library_Id" = "Library_ID")) %>%
    left_join(Probe_Sets_Core_Unit, by = c("Probe_Set" = "Name")) %>%
    mutate(
      Sample_type_submitted = "Library",
      Account_name = ifelse(is.na(Cost_Center), library.Projects, as.character(Cost_Center)),
      Contact_person_email = contact,
      Core_Unit_MPI_EVA_ID = library.CoreDB_Id,
      Customer_PandoraID = library.Full_Library_Id,
      Sampling_Protocol = "",
      Comments = "",
      Weight_of_total_Material_Provided = "",
      Weight_of_Material_to_be_used = "",
      Comments_Samples = "",
      Pre_Treatment = "",
      fluidX_Rack_ID_of_Lysate = "",
      fluidX_tube_ID_of_Lysate = "",
      Volume_of_total_Lysate_Provided = "",
      Volume_of_Lysate_to_be_used = "",
      Comments_Extraction = "",
      Pandora_batch_ID = "",
      fluidX_Rack_ID_of_Extract = "",
      fluidX_tube_ID_of_Extract = "",
      Volume_of_total_Extract_Provided = "",
      Volume_DNA_Extract_to_be_used = "",
      UDG_Treatment = "",
      Comments_LibraryPrep = "",
      Shotgun_sequencing = "no",
      Plate_ID_Pandora_batch_ID_ShotgunSequencing = "",
      Plate_position_ShotgunSequencing = "",
      Sequencing_depth_ShotgunSequencing = "",
      Sequencing_package_ShotgunSequencing = "",
      Comments_ShotgunSequencing = "",
      Capture = "yes",
      Plate_ID_Pandora_batch_ID_Capture = library.Batch,
      Plate_position_Capture = library.Position_on_Plate,
      ProbeSet = MPI_Probe_Pool, 
      Comments_Capture = "",
      Capture_Sequencing = "yes",
      Sequencing_depth_Capture = Sequencing_Depth,	
      Sequencing_package_Capture = Sequencing_Kit,	
      Comments_CaptureSequencing = "",
      Shipping = "",
      Plate_ID_Pandora_batch_ID_Aliquoting = "",
      Plate_position_Aliquoting = "",
      Destination_address = "",
      Comments_Aliquoting = ""
    ) %>%
    select(!contains("library.")) %>%
    select(!contains(c("Sequencing_Kit","Cost_Center","Lab", "MPI_Probe_Pool"))) %>%
    select(-Sequencing_Depth) %>%
    select(-Probe_Set)
  nameOutput <- paste(paste("CaptureRequest","CoreUnit",requestDate, sep = "_"),"tsv",sep = ".")
  outputTSV <- paste(opts$outdir,nameOutput, sep = "/")
  write_tsv(CoreUnitTSV, file = outputTSV)
}


