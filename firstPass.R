# Conversion to psychosis in NAPLS
# Author: A.M. Chekroud, October 2016


#### Housekeeping
# Load libraries
libs <- c("RColorBrewer", "ape", "dynamicTreeCut", "cluster",
          "doMC", "readxl", "stringr","foreach", "tidyverse")
invisible(lapply(libs, require, character.only = TRUE))
registerDoMC(detectCores()-1)

# Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("~/Documents/PhD/Projects/conversion")
workDir  <- getwd()
dataDir  <- file.path(workDir, "data") 

source("cannon_cluster_functions.R")
# list.files(dataDir)

keepVisits <- c(1:5)


pop_raw <- file.path(dataDir, "napls_pops.xlsx") %>% read_excel()
pop     <- pop_raw %>%
    filter(SubjectType == "Prodromal") %>%
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(ID, pops_criteria_met)
    


gaf_raw <- file.path(dataDir, "napls_gaf.xlsx") %>% read_excel()
gaf     <- gaf_raw %>% 
    filter(SubjectType == "Prodromal") %>%
    filter(VisitNumber %in% keepVisits) %>%
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(RID, ID, VisitNumber, gaf_CurrentScore)


# cog_raw <- file.path(dataDir, "napls_matrics.xlsx") %>% read_excel()
# cog     <- cog_raw %>% 
#     filter(SubjectType == "Prodromal") %>%
#     filter(VisitNumber %in% keepVisits) %>% 
#     mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
#     mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
#     select(RID,ID, VisitNumber, matrics_tmt_ptile)


sop_raw   <- file.path(dataDir, "napls_sops.xlsx") %>% read_excel()
sop       <- sop_raw %>% 
    filter(SubjectType == "Prodromal") %>%
    filter(VisitNumber %in% keepVisits) %>% 
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(RID, ID, VisitNumber, ends_with("_SOPS")) %>%
    mutate(Px_tot = (P1_SOPS + P2_SOPS + P3_SOPS + P4_SOPS + P5_SOPS),
           Nx_tot = (N1_SOPS + N2_SOPS + N3_SOPS + N4_SOPS + N5_SOPS + N6_SOPS),
           Gx_tot = (G1_SOPS + G2_SOPS + G3_SOPS + G4_SOPS),
           Dx_tot = (D1_SOPS + D2_SOPS + D3_SOPS + D4_SOPS)) %>%
    select(RID,ID, VisitNumber, Px_tot, Nx_tot, Gx_tot, Dx_tot)

# sop$SOP_t <- sop %>% select(ends_with("_SOPS")) %>% rowSums()
# sop       <- sop %>% select(-ends_with("_SOPS"))

gfs_raw   <- file.path(dataDir, "napls_gfs.xlsx") %>% read_excel()
gfs       <- gfs_raw %>%
    filter(SubjectType == "Prodromal") %>%
    filter(VisitNumber %in% keepVisits) %>% 
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(RID,ID, VisitNumber, gfs_current)


gfr_raw   <- file.path(dataDir, "napls_gfr.xlsx") %>% read_excel()
gfr       <- gfr_raw %>%
    filter(SubjectType == "Prodromal") %>%
    filter(VisitNumber %in% keepVisits) %>% 
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    mutate(RID = paste(SiteNumber,SubjectNumber,VisitNumber, sep="_")) %>%
    select(RID,ID, VisitNumber, gfr_current)




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


dev.off()
setEPS(); postscript("gaf_traj_by_cluster.eps")
gaf_raw %>% 
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    select(ID, VisitNumber, gaf_CurrentScore) %>%
    left_join(., select(temp, ID, cluster), by = "ID") %>% 
    mutate(visit = plyr::mapvalues(VisitNumber, 
                             from = c(1,2,3,4,5,6),
                             to = c(0,6,12,18,24,30) ) ) %>%
    group_by(visit, cluster) %>%
    summarise(meanGAF = mean(gaf_CurrentScore, na.rm=TRUE), count = n()) %>%
    ggplot() +
    geom_line(aes(x=visit, y= meanGAF, color = factor(cluster)))
dev.off()


dev.off()
png("last_gaf_by_cluster.png")
gaf_raw %>%
    mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
    select(ID, VisitNumber, gaf_CurrentScore) %>%
    left_join(., select(temp, ID, cluster), by = "ID") %>%
    arrange(ID, VisitNumber) %>%
    filter(complete.cases(.)) %>%
    group_by(ID) %>%
    summarise_each(funs(last)) %>%
    ggplot() +
    geom_violin(aes(x= factor(cluster), y = gaf_CurrentScore, fill = factor(cluster)),
                alpha=0.5) +
    geom_boxplot(aes(x= factor(cluster), y = gaf_CurrentScore),
                 width=0.1, fill="white", outlier.size=0, alpha=0.5) +
    ggtitle("Last available GAF_current score by cluster")
dev.off()


# cog_raw %>%
#     mutate(ID = paste(SiteNumber,SubjectNumber, sep="_")) %>%
#     select(ID, VisitNumber, matrics_cptip_tscore) %>%
#     left_join(., select(temp, ID, cluster), by = "ID") %>%
#     arrange(ID, VisitNumber) %>%
#     filter(complete.cases(.)) %>%
#     group_by(ID) %>%
#     summarise_each(funs(last)) %>%
#     ggplot() +
#     geom_violin(aes(x= factor(cluster), y = matrics_cptip_tscore))


# ditch matrics
## done

# gfs means gf_social gfr means gf-role

# only keep prodromes not controls
## done

# for function use gfs_current and gfr_current for each time pt
#      check which one of these has the highest range withing subject over time

## done, GFR 


# for symptoms use SOPS sum(P1-P5) (may as well also calculate N1-Nx and G1-Gx and D1-Dx but dont use them yet)
# try and include timepoints 2 and 4 now that we excluded matrics
# lets flip SOPS because currently, high SOPS means bad. So if we flip it then high SOPS will be good (consistent with the other scales). so a 0 SOPS should become a 6. ie. (SOPS * -1) + 6 
# also jam them into a 0-6 scale i.e. for gfs and gfr you can *6 and divide by 10










