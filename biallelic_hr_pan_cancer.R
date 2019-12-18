rm(list=ls()) 
library(data.table)
library(plyr)

##################################################################
load_cached  <- 0 # Switch to 1 to load from cached assembled data

##################################################################
## Gobal variables for the whole script
gene_list <-scan(file="input/Master_List.txt", what="character"); 

nGenes <- length(gene_list);
pan_cancer_types <- read.delim("input/Cancers.txt",header=F)
pan_cancer_types <- gsub(" ","",as.character(pan_cancer_types$V1))
pan_cancer_types # Output as debugging

##COAD and READ are merged under 'COLO'
pan_cancer_types <- pan_cancer_types[-c(5,24)] # Remove 2 cancer types
pan_cancer_types <- c(pan_cancer_types,"COLO")

##Read in Final_List of samples
samples_final <- read.delim("input/Final_List.txt",header=F)
samples_final <- as.character(samples_final$V1)
length(samples_final)

###################################################################
## Function to read/assemble data from raw files or cached data
load_assemble_genomics  <- function (  cancer_list, # list to loop through 
                                       path_fun, # function to generate path 
                                       cached=0, # if cached, load saved data 
                                       cached_data, # path for cached assembled data
                                       sample_select = samples_final,
                                       gene_select = gene_list,
                                       sample_ID = "V1", # sample ID from input
                                       gene_ID = "V2", # gene ID from import
                                       header=FALSE, # Whether raw data have header
                                       process_fun = function(x) x # processing function data
                                          ) {
    # Loop through cancer list and return assembled data
    if (cached) {
    load(cached_data)
    combined_data <- combined_data
    } else {
    combined_data <- list()
    for (k in cancer_list) {
      writeLines(paste("-", k))
      path <- path_fun(k)
      d <- fread(path, header=header, colClasses=c("character"))
      d <- process_fun(d) # Some magic here such as melt the data shape
      if (!is.na(gene_ID)) {
        d <- d[(d[[sample_ID]] %in% sample_select) & (d[[gene_ID]] %in% gene_select) ]
        combined_data <- rbind(d, combined_data)
      } else { # Some input data may not have gene ID
        d <- d[(d[[sample_ID]] %in% sample_select) ]
        combined_data <- rbind(d, combined_data)
      }
    }
    save(combined_data, file=cached_data)
  }
  return(combined_data)
}

  writeLines("Loading somaticMutsnps - for VUS coding") 
  somaticMutsnps <- load_assemble_genomics(pan_cancer_types, 
                  function(x) sprintf("stddata__2016_01_28/%s/20160128/Somatic_TCGA_VUS.txt", x) ,
                  cached=load_cached, cached_data="./cache/assembled_somaticMutsnps.RData",
                    process_fun=function(x) x[, .(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V43 )] ) 

  writeLines("Loading somaticMut - for LoF coding") 
  somaticMut <- load_assemble_genomics(pan_cancer_types, 
                    function(x) sprintf("stddata__2016_01_28/%s/20160128/LoF_Somatic.dat", x), 
                    cached=load_cached, cached_data="./cache/assembled_somaticMut.RData",
                     process_fun=function(x) x[, .(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V43 )] ) 

  writeLines("Loading germlineMut")
  germlineMut <- load_assemble_genomics(pan_cancer_types, 
                    function(x) sprintf("input/Germline_files/LoF_germline_Ruomu_Git_%s_Id_Gene_R.txt",x),
                    cached=load_cached, cached_data="./cache/assembled_germlineMut.RData")

  writeLines("Loading lst")
  lst  <- load_assemble_genomics(pan_cancer_types, 
              function(x) sprintf("input/LSTs/LST-%s.dat", x),
              cached=load_cached, cached_data="./cache/assembled_LST.RData",
              sample_ID="V2", gene_ID=NA )

  writeLines("Loading lohMut")
  lohMut <- load_assemble_genomics(pan_cancer_types,
              function(x) sprintf("input/LOH_calls/loh_%s_header.dat", x),
              cached=0, cached_data="./cache/assembled_lohMut.RData",
              header=TRUE,
              process_fun= function(x) melt(x, id.vars="TCGA_ID", variable.factor=FALSE), 
              sample_ID="TCGA_ID", gene_ID="variable" )

  writeLines("Loading cnvMut") # Not used in the final hrd table
  cnvMut <- load_assemble_genomics(pan_cancer_types,
              function(x) sprintf("input/CN_calls/totalcn_%s_header.dat", x),
              cached=0, cached_data="./cache/assembled_cnvMut.RData",
              header=TRUE,
              process_fun= function(x) melt(x, id.vars="TCGA_ID", variable.factor=FALSE), 
              sample_ID="TCGA_ID", gene_ID="variable" )


## Start to code with logics
hrd <- expand.grid(sample=samples_final, gene=gene_list)
hrd <- data.table(hrd)
setkey(hrd, sample, gene)


hrd$germline <- as.numeric(0) # Use default 0 as indicating no germline mut
  germlineMut$val <- 1.0
  setDT(germlineMut)
    germlineMut$sample <- germlineMut$V1
    germlineMut$gene <- germlineMut$V2
    setkey(germlineMut, sample, gene)
    germlineMut = unique(germlineMut)     # Only keep the unique ones

  hrd[germlineMut, germline := i.val]

hrd$loh <- as.double(NA)
  lohMut$value <- as.double(lohMut$value)
  setDT(lohMut)
    lohMut$sample <- lohMut$TCGA_ID
    lohMut$gene <- lohMut$variable
    setkey(lohMut, sample, gene)
    
  hrd[lohMut, loh := i.value ]

# Load and code somaticMut info
hrd$somatic <- as.double(0)
  somaticMut$sample <- somaticMut$V1
  somaticMut$gene <- somaticMut$V2
  setkey(somaticMut, sample, gene)
  somaticMut <- unique(somaticMut)
  hrd[somaticMut, somatic:=1]

########################################################
# This section deals with the mc3 somatic mutaitons
clin_sig_cat = fread("input/list_clin_sig_fields_mc3_public.csv") 
clin_sig_cat[ designated_category =="", designated_category:="vus"]

# Create the list to categorize the clin_sig type
clin_sig_pathogenic_cat = unique(clin_sig_cat[designated_category=="pathogenic", V1])
clin_sig_vus_cat = unique(clin_sig_cat[designated_category=="vus", V1])
clin_sig_benign_cat = unique(clin_sig_cat[designated_category=="benign", V1])

# Load new somatic mutation from  mc3.v0.2.8.PUBLIC.maf
if (load_cached) {
  mc3_somaticMut <- fread("cache/mc3_somatic_mutation.csv")
  mc3_somatic_vus <- fread("cache/mc3_somatic_vus.csv")
  # Previously saved cache contains other genes mutation as well
} else {
  mc3 <- fread("input/mc3.v0.2.8.PUBLIC.maf")
  mc3$sample <- substr(mc3$Tumor_Sample_Barcode, 1, 12)
  mc3 <- mc3[ sample %in% samples_final, ]
  mc3$gene <- mc3$Hugo_Symbol

  # all the variant classification cases are
  # > unique(mc3$Variant_Classification)
  #  [1] "Missense_Mutation"      "Silent"                 "5'Flank"               
  #  [4] "3'UTR"                  "RNA"                    "In_Frame_Del"          
  #  [7] "Nonsense_Mutation"      "Splice_Site"            "Intron"                
  # [10] "5'UTR"                  "In_Frame_Ins"           "Frame_Shift_Del"       
  # [13] "Nonstop_Mutation"       "3'Flank"                "Frame_Shift_Ins"       
  # [16] "Translation_Start_Site"

  mc3.LoF.list <- c( "Frame_Shift_Del",
                     "Frame_Shift_Ins",
                     "In_Frame_Del", 
                     "In_Frame_Ins",
                     "Nonsense_Mutation",
                     "Nonstop_Mutation",
                     "Splice_Site",
                     "Translation_Start_Site") 

  ## mc3 somatic lof mutation categorization
  # For mc3 somatic data, also count the pathogenic mutation annotated in the CLIN_SIG field
  mc3_somaticMut <- mc3[ (Variant_Classification %in% mc3.LoF.list) | 
                         (CLIN_SIG %in% clin_sig_pathogenic_cat) ,
                         .(sample, gene, Variant_Classification)]
  mc3_somaticMut <- mc3_somaticMut[ gene %in% gene_list,]
  fwrite(mc3_somaticMut, "cache/mc3_somatic_mutation.csv")

  ## mc3 somatic vus categorization
  mc3_somatic_vus = mc3[ !( Variant_Classification %in% mc3.LoF.list ) & 
                        !( CLIN_SIG %in% clin_sig_pathogenic_cat ) &
                        !( CLIN_SIG %in% clin_sig_benign_cat ) & 
                        ( Variant_Classification != "Silent") ,
                .(sample, gene, Variant_Classification)]
  mc3_somatic_vus <- mc3_somatic_vus[ gene %in% gene_list,]
  fwrite(mc3_somatic_vus, "cache/mc3_somatic_vus.csv")
}

  # Coding the somatic lof mutaiton in hrd object
  setkey(mc3_somaticMut, sample, gene)
  mc3_somaticMut <- unique(mc3_somaticMut)
  hrd$mc3_somatic <- 0
  hrd[mc3_somaticMut, mc3_somatic:=1]

  # coding the somatic vus in hrd object
  setkey(mc3_somatic_vus, sample, gene)
  mc3_somatic_vus <- unique(mc3_somatic_vus)
  hrd$mc3_vus <- 0
  hrd[mc3_somatic_vus, mc3_vus:=1]

# -- End of the section on MC3 data
##############################################################################

# Load and code  vus info
hrd$vus <- as.double(0)
  somaticMutsnps$sample <- somaticMutsnps$V1
  somaticMutsnps$gene <- somaticMutsnps$V2
  setkey(somaticMutsnps, sample, gene) 
  somaticMutsnps <- unique(somaticMutsnps)
  hrd[somaticMutsnps, vus:=1]

## Reading in homdel of BRCA1/2  - This is in addition to the nature comm paper
homdels <- fread("input/homdels_brca.csv")
homdels$sample <- homdels$ids 
setkey(homdels, sample, gene)
hrd[homdels, germline:=1, ] 

# Assign category per gene, per sample
hrd[, code:=.(germline*1000+ loh*100+ somatic*10 + vus )]
hrd$cat <- mapvalues(hrd$code, 
                     from = c(NA, 0, 1, 10,  11, 100, 101, 
                             110, 111, 1000,  1001, 1011, 
                             1010, # This category of 1010, could be compound pathogenic mutation
                             1101, 1100, 1111, 1110 ), 
                      to = c(NA, 'Wild_type', 'Mono_allelic_VUS', 'Mono_allelic_path', 'Mono_allelic_path', 'Wild_type', 'Biallelic_VUS',
                            'Biallelic_path', 'Biallelic_path', 'Mono_allelic_path', 'Mono_allelic_path', 'Biallelic_path',
                            'Biallelic_path', # This category of 1010, could be compound pathogenic mutation
                            'Biallelic_path', 'Biallelic_path', 'Biallelic_path', 'Biallelic_path') )

## Generate the field with mc3 somatic mutation data
hrd[, mc3_code:=.(germline*1000+ loh*100+ mc3_somatic*10 + mc3_vus )]
hrd$mc3_cat <- mapvalues(hrd$mc3_code, 
                     from = c(NA, 0, 1, 10,  11, 100, 101, 
                             110, 111, 1000,  1001, 1011, 
                             1010, # This category of 1010, could be compound pathogenic mutation
                             1101, 1100, 1111, 1110 ), 
                      to = c(NA, 'Wild_type', 'Mono_allelic_VUS', 'Mono_allelic_path', 'Mono_allelic_path', 'Wild_type', 'Biallelic_VUS',
                            'Biallelic_path', 'Biallelic_path', 'Mono_allelic_path', 'Mono_allelic_path', 'Biallelic_path',
                            'Biallelic_path', # This category of 1010, could be compound pathogenic mutation
                            'Biallelic_path', 'Biallelic_path', 'Biallelic_path', 'Biallelic_path') )

# -- End of the logic to generate hrd matrix
#######################################################################

####################################################################
#                     _
#    ____ ___  ____ _(_)___
#   / __ `__ \/ __ `/ / __ \
#  / / / / / / /_/ / / / / /
# /_/ /_/ /_/\__,_/_/_/ /_/
#
#
####################################################################

## The main logic of the script, similar to python script
main <- function() {
  ## Save the hrd object 
  fwrite(hrd, file="hrd_mc3.csv", na="NA")
}

main()
