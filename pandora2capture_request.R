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

validate_probe_set <- function(option, opt_str, value, parser) {
  valid_entries <- c("1240k","16S-rRNA","Archaic Admixture","Big Yoruba","Borrelia/Borreliella","Brucella","Candida albicans",
                     "Canid","Capra (goat)","Cowpox","Fermentation (bacteria)","Fermentation (yeast)","Helicobacter pylori",
                     "Hepatitis B","Human Screening","Immunocapture","Mitochondrial","Mitochondrial Bison",
                     "Mitochondrial_cow (Tübingen)","Mitochondrial_deer (Tübingen)","Mix Mitochondrial/Human Screening",
                     "Mycobacterium leprae","Mycobacterium tuberculosis","oralStrep-ABP","oralStrep-genomes","Ovis (sheep)",
                     "Parvovirus B19","Plasmodium Mitogenome","Plasmodium Nuclear","Plasmodium Plastid (Apicoplast)","Salmonella enterica",
                     "Staphylococcus aureus","Streptococcus mutans","Streptococcus pyogenes","Treponema pallidum",
                     "Treponema succinifaciens","Y-Chromosom","Yersinia pestis") 
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[pandora2capture_request.R] error: Invalid probe set: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

validate_lab <- function(option, opt_str, value, parser) {
  valid_entries <- c("CoreUnit","JenaLab") 
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[pandora2capture_request.R] error: Invalid lab: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

validate_kit <- function(option, opt_str, value, parser) {
  valid_entries <- c("HiSeq4000 75 SR","HiSeq4000 75 PE**","HiSeq4000 150 PE**",
                     "Flex 75 PE*", "Flex 75 SR*",
                     "MiSeqNano 150 PE", 
                     "MiSeqMicro 150 PE","MiSeq 75 PE", "MiSeq 150 PE","MiSeq 300 PE", 
                     "NextSeq500MO 75 PE", "NextSeq500MO 75 SR", 
                     "NextSeq500MO 150 SR", "NextSeq500MO 150 PE**", 
                     "NextSeq500HO 75 PE**", "NextSeq500HO 75 SR**",
                     "NextSeq500HO 150 PE**") 
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[pandora2capture_request.R] error: Invalid sequencing kit: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

validate_depth <- function(option, opt_str, value, parser) { 
  valid_entries <- c("1M","2M","5M","10M","20M", "40M")
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\n[pandora2capture_request.R] error: Invalid depth: '", value, 
                                                       "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
  }

## Parse arguments ----------------------------
parser <- OptionParser(usage = "%prog [options] .credentials")
parser <- add_option(parser, c("-f", "--file_libraries_id"), type = 'character', 
                     action = "store", dest = "file_libraries_id", 
                     help = "File containing all the libraries ID for which a capture request sheet must be produced",
                     default = NA)
parser <- add_option(parser, c("-p", "--probe_set"), type = 'character',
                     action = "callback", dest = "probe_set",
                     callback = validate_probe_set, default=NA,
                     help = "The desired probe set for the capture. To see options run: pandora2capture_request.R  -p .credentials")
parser <- add_option(parser, c("-l", "--lab"), type = 'character',
                     action = "callback", dest = "lab",
                     callback = validate_lab, default=NA,
                     help = "The lab for which the capture request is for. Should be one of: 'CoreUnit';'JenaLab'."
                     )
parser <- add_option(parser, c("-r", "--researcher"), type = 'character',
                     action = "store", dest = "researcher", default=NA,
                     help = "Name of the contact person for the request"
                     )
parser <- add_option(parser, c("-e", "--emailresearcher"), type = 'character',
                     action = "store", dest = "email", default=NA,
                     help = "Email of the contact person for the request"
                     )
parser <- add_option(parser, c("-k", "--kit_sequencing"), type = 'character',
                     action = "callback", dest = "kit",
                     callback = validate_kit, default=NA,
                     help = "The sequencing kit you will like to use. To see options run: pandora2capture_request.R  -k .credentials"
                     )
parser <- add_option(parser, c("-d", "--depth_sequencing"), type = 'character',
                     action = "callback", dest = "depth", 
                     callback = validate_depth, default="20M",
                     help = "The number of reads that you will like your captured library to be sequenced. By default the values is 20M reads. To see options run: pandora2capture_request.R  -d .credentials"
)
parser <- add_option(parser, c("-c", "--cost_center"), type = 'character',
                     action = "store", dest = "costCenter", default=NA,
                     help = "Specify the cost center. By default, this will correpond to the project in the library tab in pandora" #Default the project in Pandora, specify otherwise
)
parser <- add_option(parser, c("-o", "--outDir"), type = 'character',
                     action = "store", dest = "outdir",
                     help= "The desired output directory. By default, it is the current directory.",
                     default = "."
                     )

arguments <- parse_args(parser, positional_arguments = 1)
opts <- arguments$options

cred_file <- arguments$args
libraries <- read.csv(opts$file_libraries_id, header = F)
probe_set <- opts$probe_set
lab <- opts$lab
contact <- paste(opts$researcher, opts$email, sep = ",")
kit <- opts$kit
depth <- opts$depth
costCenter <- opts$costCenter
requestDate <- Sys.Date()
nameOutput <- paste(paste("CaptureRequest",probe_set,lab,requestDate, sep = "_"),"tsv",sep = ".")
outputTSV <- paste(opts$outdir,nameOutput, sep = "/")

con <- sidora.core::get_pandora_connection(cred_file)



libraryTab <- get_df("TAB_Library", con) %>%
  convert_all_ids_to_values(con) %>%
  filter(library.Full_Library_Id %in% libraries$V1)

if (lab == "JenaLab") {
  JenaTSV <- libraryTab %>%
    mutate(
      Library_ID = library.Full_Library_Id, 
      Researchers_Name = contact,
      Batch = library.Batch,
      Robot = case_when(
        library.Robot == FALSE ~ "n",
        library.Robot == TRUE ~ "y",
        is.na(library.Robot) ~ "NA"
      ),
      Protocol = library.Protocol,
      Probe_Set = paste0(probe_set),
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
      additional_notes = paste("Please sequence",depth,"reads in the",kit, sep = " "),
           ) %>%
    select(!contains("library."))
  write_tsv(JenaTSV, file = outputTSV)
} else if (lab == "CoreUnit") {
  Probe_Sets_Core_Unit <- read.csv("~/Documents/Postdoc/Capture_requests/Probe_Set.csv")
  code_probe_set <- Probe_Sets_Core_Unit %>%
    filter(Name == probe_set)
  CoreUnitTSV <- libraryTab %>%
    mutate(
      Sample_type_submitted = "Library",
      Account_name = case_when(
        !is.na(costCenter) ~ costCenter, 
        is.na(CostCenter) ~ library.Projects),
      Contact_person_email = contact,
      Core_Unit_MPI_EVA_ID = library.CoreDB_Id,
      Customer_PandoraID = library.Full_Library_Id,
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
      Probe_Set = code_probe_set$MPI_Probe_Pool, 
      Comments_Capture = "",
      Capture_Sequencing = "yes",
      Sequencing_depth_Capture = depth,	
      Sequencing_package_Capture = kit,	
      Comments_CaptureSequencing = "",
      Plate_ID_Pandora_batch_ID_Aliquoting = "",
      Plate_position_Aliquoting = "",
      Destination_address = "",
      Comments_Aliquoting = ""
    ) %>%
    select(!contains("library."))
  write_tsv(CoreUnitTSV, file = outputTSV)
}


