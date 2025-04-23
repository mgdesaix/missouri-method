	############################################################################################################
### Run gscramble simulations for complex 5 pops
### February 3, 2025
### Matt DeSaix
############################################################################################################

myargs <- commandArgs(trailingOnly=TRUE)
# gsp_id <- myargs[1]
start_iteration <- as.numeric(myargs[1])
final_iteration <- as.numeric(myargs[2])
outdir <- myargs[3]
ref_genos <- myargs[4]
summary_outdir <- myargs[5]

library(gscramble)
library(tidyverse)

# Load data
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
  ### Part 1a: gscramble MO pops 1,2,7
  ###
  ##############################################################################################################

  ped <- paste0(ref_genos, ".ped")
  map <- paste0(ref_genos, ".map")
  ref_gscramble <- plink2gscramble(ped = ped, map = map)
  rownames(ref_gscramble$Geno) <- ref_gscramble$I_meta$indiv
  
################## 1.1 Specify pops in GSP #############################

# Specify the not main pop (i.e. not popB)
  RepPop <- tibble(
    index = as.integer(c(1,1, 2,2)),
    pop = c("p1", "p2", "p1", "p2"),
    group = c("Pop1", "Pop7", 
              "Pop3", "Pop14")
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

ind1_options <- c("h-1-1-3-1-1", "h-1-1-3-2-2")
ind2_options <- c("h-1-2-3-1-3", "h-1-2-3-2-4")

ind1 <- sample(ind1_options, 1)
ind2 <- sample(ind2_options, 1)
ind1_index <- which(Markers$ret_ids$indiv %in% ind1)
ind2_index <- which(Markers$ret_ids$indiv %in% ind2)

ind3 <- Markers$ret_ids %>%
  filter(group == "Pop2") %>%
  slice_sample(n = 1) %>%
  pull(indiv)
ind3_index <- which(Markers$ret_ids$indiv %in% ind3)

new_index <- rep(c(ind1_index, ind2_index, ind3_index), each = 2)
ref_gscramble$Geno <- Markers$ret_geno[new_index, ]
ref_gscramble$I_meta <- Markers$ret_ids[new_index, ]
ref_gscramble$I_meta$group <- c("F1_1", "F1_1", "F1_2", "F1_2", "Pop2", "Pop2")
rownames(ref_gscramble$Geno) <- ref_gscramble$I_meta$indiv

  ##############################################################################################################
  ###
  ### Part 1b: Pedigree of the 2 F1s' offspring with a Pop2 cross
  ###
  ##############################################################################################################


  ################## 1.1 Specify pops in GSP #############################

  gsp_f1_f1b <- create_GSP("p1", "p2", F1=TRUE, F1B = TRUE)
  gsp_f1_f1b$hpop1 <- c("p1", "p2", "p3", NA, NA)
  gsp_f1_f1b$hpop2 <- c("p1", "p2", "p3", NA, NA)
  
  RepPop <- tibble(
    index = as.integer(c(1,1,1)),
    pop = c("p1", "p2", "p3"),
    group = c("F1_1", "F1_2", "Pop2")
  )

  Input_tibble <- tibble(
    gpp = list(gsp_f1_f1b),
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
    preserve_individuals = TRUE,
    preserve_haplotypes = TRUE
  )

  ################## 1.3 extract the data #############################

  # Pull all simulated individuals
  outname <- paste0(outdir, "/complex_5pops.", iteration)
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

  SysCommand <- paste0("cd ", outdir, "; admixture -j8 --supervised ", merged, ".bed ", k)
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
    add_column("iteration" = iteration,
               "pops" = "5pops") %>%
    mutate(id = as.numeric(str_split_i(Individual, "-", 4)),
           pedigree = ifelse(id == 4, "F2",
                                    ifelse(id == 5, "F2B1", NA))) %>%
    cbind(Q_mat) %>% select(-c(id))

  # Remove intermediary files
  SysCommand <- paste0("rm ", outdir, "/complex_5pops.", iteration, ".*")
  system(SysCommand)

  summary_outname <- paste0(outdir, "/complex_5pops.", iteration, ".csv")
  write_csv(Q_meta, file = summary_outname)

} # End full iteration from Part 0

file_list <- paste0(outdir, "/complex_5pops.", start_iteration:final_iteration, ".csv")
file.df <- file_list %>% lapply(read_csv) %>% bind_rows()
summary_outname <- paste0(summary_outdir, "/complex_5pops.", start_iteration, "-", final_iteration, ".csv")
write_csv(file.df, file = summary_outname)

