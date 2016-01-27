node_pos <- function(n,max.val) {
  xpos <- round(max.val*cos(2*pi*(1:n/n)),4)
  ypos <- round(max.val*sin(2*pi*(1:n/n)),4)
  results <- data.frame(id=1:n,x=xpos,y=ypos)
  return(results)
}

poly_pos <- function(scores,nodes,max.val) {
  results <- data.frame(id=numeric(),point=numeric(),x=numeric(),y=numeric())
  for(s in 1:nrow(scores)) {
    scored.nodes <- nodes[scores[s,]>0,]
    if(nrow(scored.nodes) == 1) {
      results <- rbind(results,c(s,1,scored.nodes$x,scored.nodes$y))
    } else {
      pairs <- data.frame(p1=scored.nodes$id,p2=c(scored.nodes$id[2:nrow(scored.nodes)],scored.nodes$id[1]))
      for(np in 1:nrow(pairs)) {
        score1 <- scores[s,pairs$p1[np]]
        score2 <- scores[s,pairs$p2[np]]
        score.sum <- score1 + score2
        pos.x <- ((score1*scored.nodes$x[scored.nodes$id==pairs$p1[np]])+(score2*scored.nodes$x[scored.nodes$id==pairs$p2[np]]))/score.sum
        pos.y <- ((score1*scored.nodes$y[scored.nodes$id==pairs$p1[np]])+(score2*scored.nodes$y[scored.nodes$id==pairs$p2[np]]))/score.sum
        results <- rbind(results,c(s,np,pos.x,pos.y))
      }
    }
  }
  
  names(results) <- c("id","point","x","y")
  results <- results[order(results$id),]
  return(results)
}

centroid_pos <- function(polys) {
  results <- data.frame(id=numeric(),x=numeric(),y=numeric())
  for(i in unique(polys$id)) {
    poly <- polys[polys$id == i,]
    if(nrow(poly) == 1) {
      c.x <- poly$x
      c.y <- poly$y
    }
    if(nrow(poly) == 2) {
      c.x <- (min(poly$x)+max(poly$x))/2
      c.y <- (min(poly$y)+max(poly$y))/2
    }
    if(nrow(poly) > 2) {
      poly <- rbind(poly,poly[1,])
      a.sum <- 0
      x.sum <- 0
      y.sum <- 0
      for(p in 1:(nrow(poly)-1)) {
        a.sum <- a.sum + (poly$x[p] * poly$y[p + 1] - poly$x[p+1] * poly$y[p])
        x.sum <- x.sum + (poly$x[p] + poly$x[p + 1]) * (poly$x[p] * poly$y[p + 1] - poly$x[p + 1] * poly$y[p])
        y.sum <- y.sum + (poly$y[p] + poly$y[p + 1]) * (poly$x[p] * poly$y[p + 1] - poly$x[p + 1] * poly$y[p])
      }
      a <- 0.5*abs(a.sum)
      c.x <- x.sum/(6*a)
      c.y <- y.sum/(6*a)
    }
    results <- rbind(results,c(i,c.x,c.y))  
  }
  names(results) <- c("id","x","y")
  return(results)
}

plot_landscape <- function(votes,clusters,positions=NULL) {
  library(dplyr)
  library(ggplot2)
  
  votes_data <- votes %>%
    filter(primvec %in% clusters) %>%
    select(one_of(clusters))
  
  most_votes <- max(votes_data)
  
  nodes <- node_pos(length(clusters),most_votes)
  polygons <- poly_pos(votes_data,nodes,most_votes)
  centroids <- centroid_pos(polygons)
  plot_centroids <- centroids %>% group_by(x,y) %>% summarize(z = n())
  
  p <- ggplot(centroids,aes(x=x,y=y)) + 
    stat_density2d(geom="polygon") +
    geom_jitter(color="black",alpha=0.5,position = position_jitter(height=0.2,width=0.2),size=0.7) +
    xlim(-2*most_votes,2*most_votes) + ylim(-2*most_votes,2*most_votes) +
    theme_classic()
  
  return(p)
}
