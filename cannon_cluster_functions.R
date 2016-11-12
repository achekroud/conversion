clusteringPipeline <- function(inputMatrix, distMethod = "manhattan", hclustMethod = "ward.D2",
                               scaled = FALSE, leaf.labels, cut_height = NULL, perm.test = FALSE, nPerms = 100){
    ## This is one mega function that encompasses an entire clustering pipeline
    ## The pipeline will calculate a number of steps as follows:
    ##       Scale data (optional)
    ##       Calculate a distance matrix
    ##       Cluster the data
    ##       Cut the resultant dendrogram at the optimal points
    ##       Summarise the metrics used in this pipeline, and the results
    ##       Output the colnames for each resultant cluster
    ##       Calculate silhouette distances for the solution
    ##       Optionally conduct a permutation test using the simprof method
    ##       Return all of the above as a list
    ## Required inputs:
    ##       inputMatrix is a matrix of your to-be-clustered data
    ##       distMethod is the desired distance metric
    ##       hclustMethod is the desired linkage function
    ##       scaled is an optional flag to z-score each variable in the matrix
    ##       leaf labels can be provided as a vector
    ##         if given, they will be inserted into the dendrogram for nicer plotting
    ##       cut_height is an optional flag to insert alternative cut heights in the dendrogram  
    ##       perm.test is an optional flag to conduct a permutation test
    ##         nPerms is the number of permutations to conduct -- 100 is default.
    ##            - be careful what u wish 4, this is slooooow
    ## Outputs: a list
    
    if (scaled == TRUE) {inputMatrix = scale(inputMatrix); print("Data have been scaled")} else print("Data NOT scaled")
    print("calculating distance matrix...")
    dm  <- dist(t(inputMatrix), method = distMethod)
    print("clustering...")
    hc  <- hclust(dm, method = hclustMethod)
    cut <- cutreeDynamic(hc, minClusterSize = 1, distM = as.matrix(dm), verbose=4, cutHeight = cut_height)
    if(0 %in% cut){cut <- cut + 1}
    summary <- paste(paste0("Scaling = ", scaled),
                     paste0("Distance Metric = ", distMethod),
                     paste0("H-Clust Method = ", hclustMethod), sep = " \n  ")
    if (!missing(leaf.labels)) {hc$labels <- as.character(leaf.labels)}
    nClust <- length(table(cut))
    print(paste0(nClust, " clusters identified"))
    cNames <- list(NULL)
    for (i in 1:nClust){
        cNames [[i]] <- colnames(inputMatrix)[cut==i]
        names(cNames)[[i]] <- paste0("cluster", i)
    }
    print("calculating silhouette distances...")
    sil <- silhouette(cut, dist = dm)
    
    if (perm.test == TRUE) {
        if (!missing(leaf.labels)) {colnames(inputMatrix) <- leaf.labels}
        print("Beginning Permutation Testing...updates per 50 perms")
        perm <- simprof(inputMatrix, method.cluster = hclustMethod, method.distance = distMethod,
                        sample.orientation = "column", num.expected = nPerms, num.simulated = nPerms,
                        silent = FALSE, increment = nPerms/2)
        print("Permutation testing complete. See output$perm for more details")
        out <- list(summary = summary, dm = dm, hc = hc, cut = cut, cNames = cNames, sil = sil, perm = perm)
    } else { 
        print("No permutation testing conducted")
        out <- list(summary = summary, dm = dm, hc = hc, cut = cut, cNames = cNames, sil = sil, perm = "not_conducted")
    }
    print(summary)
    return(out)
}

plotDendro <- function(pipeline, alt.cex = 0.8, nice.pal = brewer.pal(9,"Set1")[pipeline$cut], edgeWidth = 1,
                       labelOffset = 100, title = list(pipeline$summary, cex=0.9)){
    ## This function will return a dendrogram for a given clustering object returned by clusteringPipeline
    ## Inputs:
    ##   required: pipeline object
    ##   optional: alt.cex is a scaling factor for the leaf labels
    ##             nice.pal defaults to a 9-colour set from brewer pal
    ##                 -feel free to provide your own palette if you like
    
    #### TODO : this breaks if the leaf labels dont match up with the colnames of the original matrix
    plt <- plot(ape::as.phylo(pipeline$hc), tip.color=nice.pal, label.offset=labelOffset,
                cex = alt.cex, font = 1, main = title, edge.width = edgeWidth)
}