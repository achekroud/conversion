# Conversion to psychosis in NAPLS
# Author: A.M. Chekroud, October 2016


#### Housekeeping

# Load libraries
libs <- c("RColorBrewer", "ape", "dynamicTreeCut", "cluster",
          "doMC", "readxl", "stringr","foreach", "tidyverse")
invisible(lapply(libs, require, character.only = TRUE))
registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("~/Documents/PhD/Projects/conversion")
workDir  <- getwd()
dataDir  <- file.path(workDir, "data") 

source("cannon_cluster_functions.R")
# list.files(dataDir)

keepVisits <- c(1,3,5)


pop_raw <- file.path(dataDir, "napls_pops.xlsx") %>% read_excel()
pop     <- pop_raw %>%
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(ID, pops_criteria_met)
    


gaf_raw <- file.path(dataDir, "napls_gaf.xlsx") %>% read_excel()
gaf     <- gaf_raw %>% 
    filter(VisitNumber %in% keepVisits) %>%
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(RID, ID, VisitNumber, gaf_CurrentScore)


cog_raw <- file.path(dataDir, "napls_matrics.xlsx") %>% read_excel()
cog     <- cog_raw %>% 
    filter(VisitNumber %in% keepVisits) %>% 
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(RID,ID, VisitNumber, matrics_tmt_ptile)


sop_raw   <- file.path(dataDir, "napls_sops.xlsx") %>% read_excel()
sop       <- sop_raw %>% 
    filter(VisitNumber %in% keepVisits) %>% 
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(RID, ID, VisitNumber, ends_with("_SOPS"))

sop$SOP_t <- sop %>% select(ends_with("_SOPS")) %>% rowSums()
sop       <- sop %>% select(-ends_with("_SOPS"))







all <- full_join(gaf, cog, by = "RID") %>% 
    full_join(., sop, by = "RID") %>%
    select(RID, ID, VisitNumber, gaf_CurrentScore, matrics_tmt_ptile,SOP_t) %>%
    gather(key = "measure", value = "score", gaf_CurrentScore:SOP_t) %>% 
    transmute(ID = ID, 
              score = score,
              visitmeasure = paste(VisitNumber, measure, sep="_")) %>%
    arrange(ID) %>%
    spread(key = "visitmeasure", value = "score")

all$nas <- apply(all[,2:9], 1, function(x) length(x[is.na(x)]))

all <- all %>% filter(nas < 3)

clusterMe <- as.matrix(all[,2:10])
rownames(clusterMe) <- all$ID

# corMat <- 1 - cor(t(clusterMe), use="pairwise")
# distMat <- dist(clusterMe, method = "euclidean")
# distMat <- as.dist(corMat)
# hc <- hclust(distMat, method = "ward.D")

out <- clusteringPipeline(t(clusterMe))
plotDendro(out)


str(out)

all$cluster <- out$cut

temp <- full_join(all, pop, by = "ID")

gaf_raw %>% 
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    select(ID, VisitNumber, gaf_CurrentScore) %>%
    left_join(., select(temp, ID, cluster), by = "ID") %>% 
    group_by(VisitNumber, cluster) %>%
    summarise(meanGAF = mean(gaf_CurrentScore, na.rm=TRUE)) %>%
    ggplot() +
    geom_line(aes(x=VisitNumber, y= meanGAF, color = factor(cluster)))

gaf_raw %>%
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    select(ID, VisitNumber, gaf_CurrentScore) %>%
    left_join(., select(temp, ID, cluster), by = "ID") %>%
    arrange(ID, VisitNumber) %>%
    filter(complete.cases(.)) %>%
    group_by(ID) %>%
    summarise_each(funs(last)) %>%
    ggplot() +
    geom_violin(aes(x= factor(cluster), y = gaf_CurrentScore))

cog_raw %>%
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    select(ID, VisitNumber, matrics_cptip_tscore) %>%
    left_join(., select(temp, ID, cluster), by = "ID") %>%
    arrange(ID, VisitNumber) %>%
    filter(complete.cases(.)) %>%
    group_by(ID) %>%
    summarise_each(funs(last)) %>%
    ggplot() +
    geom_violin(aes(x= factor(cluster), y = matrics_cptip_tscore))


