require(data.table)
require(dbscan)

viirs = fread('viirs_2012.csv')

clust = hdbscan(viirs[,c('x','y')], minPts = 5)

# cluster membership assignment for each
clust$cluster
# cluster membership probability of each
clust$membership_prob
# outlier score for each point, *not* complement of membership scor
clust$outlier_scores
# cluster score for each *cluster*
clust$cluster_scores

# plot hierarchical tree
plot(clust)

library(dbscan)

#dbscan
db.clust = dbscan(lung.cancer2[,c('x','y')], minPts = 5, eps=0.05)
# cluster membership assignment for each
clust$cluster
# cluster membership probability of each
clust$membership_prob
# outlier score for each point, *not* complement of membership scor
clust$outlier_scores
# cluster score for each *cluster*
clust$cluster_scores

hullplot(lung.pts,db.clust)

# hdbscan
hdb.clust = hdbscan(lung.cancer2[,c('x','y')], minPts = 5)
hullplot(lung.pts,hdb.clust)

# plot hierarchical tree
plot(hdb.clust)
plot(hdb.clust, show_flat = TRUE)

plot(lung.pts, col=clust$cluster+1L, cex = .5)


library("dbscan")
data("moons")
plot(moons, pch=20)
cl <- hdbscan(moons, minPts = 5)
cl
plot(moons, col=cl$cluster+1, pch=20)
