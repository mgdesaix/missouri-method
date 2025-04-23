############################################################################################################
### Run gscramble simulations for ghosts
### October 3, 2024
### Matt DeSaix
############################################################################################################

myargs <- commandArgs(trailingOnly=TRUE)
gsp_id <- myargs[1]
start_iteration <- as.numeric(myargs[2])
final_iteration <- as.numeric(myargs[3])
popA <- myargs[4]
outdir <- myargs[5]
ref_genos <- myargs[6]
ghost_genos <- myargs[7]
summary_outdir <- myargs[8]

library(gscramble)
library(tidyverse)


# Load data
gsp_file <- paste0("/lustrefs/nwrc/projects/DeSaix/MO-method/data/gsps/mo.", gsp_id, ".gsp.csv")
gsp <- read_csv(gsp_file, show_col_types = F)
gsp_f1 <- create_GSP("p1", "p2", F1=TRUE)


RecRates <- read_csv("/lustrefs/nwrc/projects/DeSaix/MO-method/data/gsps/MOMethod-SpatialPops1672x27866.recrates.csv", show_col_types = FALSE)
# map <- paste0(ref_genos, ".map")
# RecRates <- plink_map2rec_rates(map = map)

ref_pops <- read_table(paste0(ref_genos, ".fam"),
                       col_names = c("Pop", "Individual", "a", "b", "c", "d"),
                       show_col_types = FALSE) %>%
  pull(Pop) %>%
  unique()
k <- length(ref_pops)

##############################################################################################################
###
### Part 0: Start iteration
###
##############################################################################################################

# Number of iterations is provided by R script input
for(iteration in start_iteration:final_iteration){

  ##############################################################################################################
  ###
  ### Part 1a: gscramble MO popA
  ###
  ##############################################################################################################
  popABC <- read_delim("./data/gsps/scrambled-pop-array.txt", delim = "\t", col_names = c("A", "B", "C")) %>%
    filter(A == popA)
  popB <- popABC$B
  popC <- popABC$C

  ped <- paste0(ref_genos, ".ped")
  map <- paste0(ref_genos, ".map")
  ref_gscramble <- plink2gscramble(ped = ped, map = map)
  rownames(ref_gscramble$Geno) <- ref_gscramble$I_meta$indiv
  
################## 1.1 Specify pops in GSP #############################

# Specify the not main pop (i.e. not popA)
  RepPop <- tibble(
    index = as.integer(c(1,1)),
    pop = c("p1", "p2"),
    group = c(popB, popC)
  )

  Input_tibble <- tibble(
    gpp = list(gsp_f1),
    reppop = list(RepPop)
  )

  ################## 1.2 scramble the genotypes  #############################

  Segments <- segregate(
    request = Input_tibble,
    RR = RecRates,
    MM = ref_gscramble$M_meta
  )

  Markers <- segments2markers(
    Segs = Segments,
    Im = ref_gscramble$I_meta,
    Mm = ref_gscramble$M_meta,
    G = ref_gscramble$Geno,
    preserve_individuals = FALSE
  )

  ################## 1.3 extract the data #############################

  # Pull all simulated individuals
  outname <- paste0(outdir, "/", popA, ".ghost.", gsp_id, ".", iteration, ".scrambled")
  sim_index <- which(Markers$ret_ids$group == popA)

  gscramble2plink(I_meta = Markers$ret_ids[sim_index,],
                  M_meta = ref_gscramble$M_meta,
                  Geno = Markers$ret_geno[sim_index,],
                  prefix = outname)

  SysCommand <- paste0("plink --noweb --file ", outname, " --make-bed --out ", outname)
  system(SysCommand)

  ghost_input <- paste0(outname, ".with_ghosts")
  SysCommand <- paste0("plink --noweb --bfile ", outname, " --bmerge ", ghost_genos, " --recode --out ", ghost_input)
  system(SysCommand)
  

  ##############################################################################################################
  ###
  ### Part 1b: gscramble with ghosts
  ###
  ##############################################################################################################

  ped <- paste0(ghost_input, ".ped")
  map <- paste0(ghost_input, ".map")
  ref_gscramble <- plink2gscramble(ped = ped, map = map)
  rownames(ref_gscramble$Geno) <- ref_gscramble$I_meta$indiv


  sample_founders <- ref_gscramble$I_meta %>% group_by(group) %>% slice_sample(n=8) %>% pull(indiv)
  sample_index <- which(ref_gscramble$I_meta$indiv %in% sample_founders)

  ################## 1.1 Specify pops in GSP #############################

  RepPop <- tibble(
    index = as.integer(c(1,1)),
    pop = c("p1", "p2"),
    group = c(popA, "-")
  )

  Input_tibble <- tibble(
    gpp = list(gsp),
    reppop = list(RepPop)
  )

  ################## 1.2 scramble the genotypes  #############################

  Segments <- segregate(
    request = Input_tibble,
    RR = RecRates,
    MM = ref_gscramble$M_meta
  )

  Markers <- segments2markers(
    Segs = Segments,
    Im = ref_gscramble$I_meta[sample_index,],
    Mm = ref_gscramble$M_meta,
    G = ref_gscramble$Geno[sample_index,],
    preserve_individuals = TRUE
  )

  ################## 1.3 extract the data #############################

  # Pull all simulated individuals
  outname <- paste0(outdir, "/", popA, ".ghost.", gsp_id, ".", iteration)
  sim_index <- which(Markers$ret_ids$group == "ped_hybs")
  Markers$ret_ids$group <- "-"

  gscramble2plink(I_meta = Markers$ret_ids[sim_index,],
                  M_meta = ref_gscramble$M_meta,
                  Geno = Markers$ret_geno[sim_index,],
                  prefix = outname)

  SysCommand <- paste0("plink --noweb --file ", outname, " --make-bed --out ", outname)
  system(SysCommand)

  # merge simulated genotypes with reference
  merged <- paste0(outname, ".merged")
  SysCommand <- paste0("plink --noweb --bfile ", outname, " --bmerge ", ref_genos,
                       " --make-bed --out ", merged)
  system(SysCommand)

  ##############################################################################################################3
  ###
  ### Part 2: Admixture
  ###
  ##############################################################################################################

  # ADMIXTURE requires a *.pop file to complement *.bed/*.bim/*.fam files in supervised analyses
  SysCommand <- paste0("cut -f1-2 -d' ' ", merged, ".fam > ", merged, ".pop")
  system(SysCommand)

  SysCommand <- paste0("cd ", outdir, "; admixture -j4 --supervised ", merged, ".bed ", k)
  system(SysCommand)

  merged_fam <- read_table(paste0(merged, ".fam"),
                         col_names = c("Pop", "Individual", "a", "b", "c", "d"),
                         show_col_types = FALSE) %>%
    select(Pop, Individual)
  merged_index <- which(merged_fam$Pop == "-")

  Q_file <- paste0(merged, ".", k, ".Q")
  Q_mat <- read_table(Q_file,
                      col_names = ref_pops,
                      show_col_types = FALSE)[merged_index,]

  Q_meta <- merged_fam[merged_index,] %>%
    select(Individual) %>%
    add_column("gsp" = gsp_id,
               "iteration" = iteration,
               "popA" = popA,
               "popB" = "Ghost") %>%
    mutate(id = as.numeric(str_split_i(Individual, "-", 4)),
           pedigree = ifelse(gsp == "f1b", "F1B1", ifelse(gsp == "f1", "F1",
                             ifelse(id == 5, "F1",
                                    ifelse(id == 6, "F1B1",
                                           ifelse(id == 7, "F1B2", NA)))))) %>%
    cbind(Q_mat) %>% select(-c(gsp, id))

  # Remove intermediary files
  SysCommand <- paste0("rm ", outdir, "/", popA, ".ghost.", gsp_id, ".", iteration, ".*")
  system(SysCommand)

  summary_outname <- paste0(outdir, "/", popA, ".ghost.", gsp_id, ".", iteration, ".csv")
  write_csv(Q_meta, file = summary_outname)

} # End full iteration from Part 0

# summary_pattern <- paste0(popA, ".ghost.", gsp_id, ".*.csv")
# file_list <- list.files(path = outdir, pattern = summary_pattern, full.names = TRUE)
file_list <- paste0(outdir, "/", popA, ".ghost.", gsp_id, ".", start_iteration:final_iteration, ".csv")
file.df <- file_list %>% lapply(read_csv) %>% bind_rows()
summary_outname <- paste0(summary_outdir, "/", popA, ".ghost.", gsp_id, ".", start_iteration, "-", final_iteration, ".csv")
write_csv(file.df, file = summary_outname)

