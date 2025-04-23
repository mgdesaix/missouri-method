############################################################################################################
### Run gscramble simulations for Missouri Method: scramble reference pops
### October 3, 2024
### Matt DeSaix
############################################################################################################

myargs <- commandArgs(trailingOnly=TRUE)
popA <- myargs[1]
popB <- myargs[2]
popC <- myargs[3]
prefix <- myargs[4]
outdir <- myargs[5]
iterations <- as.numeric(myargs[6])

library(gscramble)
library(tidyverse)


# Load data
gsp <- create_GSP("p1", "p2", F1 = TRUE)
ped <- paste0(prefix, ".ped")
map <- paste0(prefix, ".map")
ref_gscramble <- plink2gscramble(ped = ped, map = map)
rownames(ref_gscramble$Geno) <- ref_gscramble$I_meta$indiv

# RecRates <- read_csv("/lustrefs/nwrc/projects/DeSaix/MO-method/data/gsps/MOMethod-SpatialPops1589.recrates.csv", show_col_types = FALSE)
RecRates <- plink_map2rec_rates(map = map)

  ##############################################################################################################3
  ###
  ### Part 1: gscramble
  ###
  ##############################################################################################################

I_meta.tmp <- as.data.frame(matrix(NA, nrow = iterations, ncol = ncol(ref_gscramble$I_meta)))
colnames(I_meta.tmp) <- colnames(ref_gscramble$I_meta)
Geno.tmp <- as.data.frame(matrix(NA, nrow = iterations, ncol = ncol(ref_gscramble$Geno)))
colnames(Geno.tmp) <- colnames(ref_gscramble$Geno)

for(i in 1:iterations){

  ################## 1.1 Specify pops in GSP #############################

# Specify the not main pop (i.e. not popA)
  RepPop <- tibble(
    index = as.integer(c(1,1)),
    pop = c("p1", "p2"),
    group = c(popB, popC)
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
    Im = ref_gscramble$I_meta,
    Mm = ref_gscramble$M_meta,
    G = ref_gscramble$Geno,
    preserve_individuals = FALSE
  )

  ################## 1.3 extract the data #############################

  # randomly sample a single individual
  sim_index <- sample(which(Markers$ret_ids$group == popA),1)
  I_meta.tmp[i,] <- Markers$ret_ids[sim_index,]
  Geno.tmp[i,] <- Markers$ret_geno[sim_index,]


} # end for loop
  outname <- paste0(outdir, "/", popA, ".scrambled.", iterations)

  gscramble2plink(I_meta = I_meta.tmp,
                  M_meta = ref_gscramble$M_meta,
                  Geno = Geno.tmp,
                  prefix = outname)

  SysCommand <- paste0("plink --noweb --file ", outname, " --make-bed --out ", outname)
  system(SysCommand)


 