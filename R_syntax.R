#Author1=ORCID:0000-0002-3324-2255
#Author2=ORCID:0000-0002-1319-0785
#Title=Effect of inaccurate edge information on centrality ranking error in classroom networks
#ver=12.221
#date=23th of May, 2025

#01.ENVIRONMENT and PACKAGES
#01.ENVIRONMENT and PACKAGES
rm(list = ls())
library("igraph")
library("readr")
library("googledrive")

#checking if kendall-tau-beta calculates correctly in R
ref1<-c(1:5)
dist1<-c(c(1,2,3,5,4)) 
cor(ref1, dist1 , use = "complete.obs", method = "kendall")
(1-cor(ref1, dist1 , use = "complete.obs", method = "kendall"))/2
#good, it calculates the way it shall
rm(ref1)
rm(dist1)

#02.LOAD THE GRAPH, RG
#02.LOAD THE GRAPH, RG

#reading the data permission given by Fire M.
#source:Fire, M., Katz, G., Elovici, Y., Shapira, B., & Rokach, L. (2012). 
#Predicting student exam’s scores by analyzing social network data. In Active Media Technology: 
#8th International Conference, AMT 2012, Macau, China, December 4-7, 2012. 
#Proceedings 8 (pp. 584-595). 
#Springer Berlin Heidelberg.

drive_download(file=as_id("https://drive.google.com/file/d/1Df4F_uIEHQAeG0fHkiyJyZZ-2jWtBx9d/")) #downloading the data
unzip(zipfile="students.zip", file="multigraph_hashAnonymized.csv") #extracting the data

cnss <- read_table("multigraph_hashAnonymized.csv", col_names = FALSE, 
                    col_types = cols(X1 = col_integer(), 
                                     X2 = col_integer(), X3 = "skip"))

colnames(cnss) <- c("node1", "node2")


#creating RG, the real, baseline graph
rg <- graph_from_data_frame(cnss, directed = FALSE)
rg <- simplify(rg, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = igraph_opt("edge.attr.comb"))
comp <- components(rg)
rg <- decompose(rg, max.comps = 1, min.vertices = max(comp$csize))[[1]] #selecting the biggest component
rm(comp)
rm(cnss)

V(rg)$name <- c(1:vcount(rg))

#saveRDS(rg, file="RGsav.RDS")

#checking the basic parameters of the real graph
vcount(rg)
ecount(rg)
edge_density(rg, loops=FALSE)
#plot(rg, layout=layout_with_fr(rg), vertex.size=5, vertex.label.cex=0.3)
table(degree(rg))
#hist(degree(rg))
median(degree(rg))
mean(degree(rg))


#ground truth centrality parameters for each node
#real graph centrality measures in a data.frame

  rgdc <- degree(rg, normalized=TRUE)
  rgrc <- (ego_size(rg, 2)-1)/(vcount(rg)-1)
  rgxc <- eccentricity(rg)
  
  rgcc <- closeness(rg)
  rgbc <- betweenness(rg)
  
  rgec <- eigen_centrality(rg)$vector
  rgkc <- authority.score(rg)$vector
  rgpc <- power_centrality(rg, rescale = TRUE , exponent = 0.9, tol = 1e-03)
  
  rghc <- harmonic_centrality(rg)
  rgpr <- page_rank(rg)$vector

#two basic parameter of the graph
rg_ec1 <- ecount(rg) #edge count of the real graph
rg_nc1 <- vcount(rg) #node count of the real graph

rg_ec1
rg_nc1

#03.SIMULATION 
#03.SIMULATION 
for (k in 1:1024)
 for (me in 0 : 52 )
  for (mte in 0 : 52 )
    {
    
    #first removing the edges
    if (me==0) {sg2<-rg} else #if edges are not removed keep the original graph, else if removed, remove
        {
         deids <- sample(x=c(1:rg_ec1), size=me, replace=FALSE) #ids of edges to delete
         sg2 <- delete_edges(rg, deids) #delete
        }
    
    #second add random edges
    if (!mte==0) {
                  for (c in 1 : mte ) #now let' add new ones
                      {
                       nodes <- sample(V(sg2), size = 2)
                       
                       while (are.connected(sg2, nodes[1], nodes[2])) 
                            { 
                            nodes <- sample(V(sg2), size = 2) 
                            } 
                       
                       sg2 <- add_edges(sg2, c(nodes[1], nodes[2]) )
                       
                      }
                  }
    
    #sg2 is the sampled, moisy, distorted graph
    #calculate centralities for sg2
    #sometimes it runs to an error, hence we use try
    #sgc = sampled graph centrality measures data frame
    
    #here we need error management, because sometimes it comes back with a memory error

      dc <- try(degree(sg2, normalized=TRUE))
      rc <- try( (ego_size(sg2, 2)-1)/(vcount(sg2)-1) ) 
      xc <- try(eccentricity(sg2))
      
      cc <- try(closeness(sg2))
      bc <- try(betweenness(sg2))
      
      ec <- try(eigen_centrality(sg2)$vector)
      kc <- try(authority.score(sg2)$vector)
      pc <- try(power_centrality(sg2, rescale = TRUE , exponent = 0.9, tol = 1e-03))
      
      hc <- try(harmonic_centrality(sg2))
      pr <- try(page_rank(sg2)$vector)
    
      while (typeof(dc)=="character") { Sys.sleep(50) ; dc <- try(degree(sg2, normalized=TRUE)) }
      while (typeof(rc)=="character") { Sys.sleep(50) ; rc <- try( (ego_size(sg2, 2)-1)/(vcount(sg2)-1) )  }  
      while (typeof(xc)=="character") { Sys.sleep(50) ; xc <- try(eccentricity(sg2)) } 
      
      while (typeof(cc)=="character") { Sys.sleep(50) ; cc <- try(closeness(sg2)) }
      while (typeof(bc)=="character") { Sys.sleep(50) ; bc <- try(betweenness(sg2)) } 
      
      while (typeof(ec)=="character") { Sys.sleep(50) ; ec <- try(eigen_centrality(sg2)$vector) }
      while (typeof(kc)=="character") { Sys.sleep(50) ; kc <- try(authority.score(sg2)$vector) } 
      while (typeof(pc)=="character") { Sys.sleep(50) ; pc <- try(power_centrality(sg2, rescale = TRUE , exponent = 0.9, tol = 1e-01)) } 
        
      while (typeof(hc)=="character") { Sys.sleep(50) ; hc <- try(harmonic_centrality(sg2)) }
      while (typeof(pr)=="character") { Sys.sleep(50) ; pr <- try(page_rank(sg2)$vector) }       
    
    #nas are replaced by zeros
    if (sum(is.na(dc))>0) {dc[is.na(dc)] <- 0}
    if (sum(is.nan(dc))>0) {dc[is.nan(dc)] <- 0}
    if (sum(is.na(rc))>0) {rc[is.na(rc)] <- 0}
    if (sum(is.nan(rc))>0) {rc[is.nan(rc)] <- 0}  
    if (sum(is.na(xc))>0) {xc[is.na(xc)] <- 0}
    if (sum(is.nan(xc))>0) {xc[is.nan(xc)] <- 0}
      
    if (sum(is.na(cc))>0) {cc[is.na(cc)] <- 0}
    if (sum(is.nan(cc))>0) {cc[is.nan(cc)] <- 0}
    if (sum(is.na(bc))>0) {bc[is.na(bc)] <- 0}
    if (sum(is.nan(bc))>0) {bc[is.nan(bc)] <- 0}  
    
    if (sum(is.na(ec))>0) {ec[is.na(ec)] <- 0}
    if (sum(is.nan(ec))>0) {ec[is.nan(ec)] <- 0}
    if (sum(is.na(kc))>0) {kc[is.na(kc)] <- 0}
    if (sum(is.nan(kc))>0) {kc[is.nan(kc)] <- 0}  
    if (sum(is.na(pc))>0) {pc[is.na(pc)] <- 0}
    if (sum(is.nan(pc))>0) {pc[is.nan(pc)] <- 0}
    
    if (sum(is.na(hc))>0) {hc[is.na(hc)] <- 0}
    if (sum(is.nan(hc))>0) {hc[is.nan(hc)] <- 0}
    if (sum(is.na(pr))>0) {pr[is.na(pr)] <- 0}
    if (sum(is.nan(pr))>0) {pr[is.nan(pr)] <- 0}  
    

    #calculating the error measure
    #first every value is 0.5 the baseline error
    tmp_dc <- 0.5
    tmp_rc <- 0.5
    tmp_xc <- 0.5
  
    tmp_cc <- 0.5  
    tmp_bc <- 0.5
    
    tmp_ec <- 0.5
    tmp_kc <- 0.5
    tmp_pc <- 0.5
    
    tmp_hc <- 0.5
    tmp_pr <- 0.5
    
    #if there is no more than 2 edge in the graph, than it is 0.5
    if ( ecount(sg2)>1 ) 
        {
      
      tmp_dc <- 0.5 *  (1 - cor(rgdc, dc , use = "complete.obs", method = "kendall"))
      tmp_rc <- 0.5 *  (1 - cor(rgrc, rc , use = "complete.obs", method = "kendall"))
      tmp_xc <- 0.5 *  (1 - cor(rgxc, xc , use = "complete.obs", method = "kendall"))
      
      tmp_cc <- 0.5 *  (1 - cor(rgcc, cc, use = "complete.obs", method = "kendall"))
      tmp_bc <- 0.5 *  (1 - cor(rgbc, bc , use = "complete.obs", method = "kendall"))
      
      tmp_ec <- 0.5 *  (1 - cor(rgec, ec, use = "complete.obs", method = "kendall"))
      tmp_kc <- 0.5 *  (1 - cor(rgkc, kc, use = "complete.obs", method = "kendall"))
      tmp_pc <- 0.5 *  (1 - cor(rgpc, pc, use = "complete.obs", method = "kendall"))
      
      tmp_hc <- 0.5 *  (1 - cor(rghc, hc, use = "complete.obs", method = "kendall"))
      tmp_pr <- 0.5 *  (1 - cor(rgpr, pr, use = "complete.obs", method = "kendall"))
        }
    
  #results are written in a files by k
  filname <- paste0("xGRAPH_20230522_k",k,".csv")

  Sys.sleep(0.01)
  
  cat(k, me, mte 
      , tmp_dc , tmp_rc , tmp_xc 
      , tmp_cc  , tmp_bc  
      , tmp_ec , tmp_kc ,  tmp_pc  
      , tmp_hc , tmp_pr , "1 \n" 
      , sep=";", file = filname, append = TRUE)
    
  Sys.sleep(0.01)
    
} 


#04. ANALYSES -- READING THE DATA
#04. ANALYSES -- READING THE DATA
rm(list = ls())
library("readr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("data.table")
library("simpleboot")


#i vector is name, is for the centrality scores
i <- c(1:10)
names(i) <- c("Degree", "Reach", "Eccentricity", "Closeness", "Betwenness",
              "Eigenvector", "Kleinberg_authority", "Power", "Harmonic", "Page_rank")

#k is the number of cycles, which we analye
k <- 1024

#edges count of RG
ec1 <- 256

res <- NULL
#reading the first file
res <- read_delim("xGRAPH_20230522_k1.csv", 
                  delim = ";", escape_double = FALSE, col_names = FALSE, 
                  col_types = cols(X1 = col_integer(), 
                                   X2 = col_integer(), X3 = col_integer(), 
                                   X14 = col_skip()), 
                  skip = ifelse(tail(readLines("GRAPH_20230520_k1.csv"), 1) == "", 1, 0) , 
                  trim_ws = TRUE)

colnames(res) <- c("k", "me", "mte"
          ,"dc" , "rc" ,"xc" 
          , "cc"  , "bc"  
          , "ec" , "kc" ,  "pc"  
          , "hc" , "pr")

res <- data.table(res)

#checking the structure of the first file
table(res$k) #it shall be a constant
table(res$me) #it shall be constant
table(res$mte) #it shall be constant
summary(res)

res1 <- NULL
res1 <- data.table()

#cycle for reading all of the rest
for (c in 2:k)
  {
  
  filname <- paste0("xGRAPH_20230522_k",k,".csv")
  
  res1 <- read_delim(file=filname, 
                    delim = ";", escape_double = FALSE, col_names = FALSE, 
                    col_types = cols(X1 = col_integer(), 
                                     X2 = col_integer(), X3 = col_integer(), 
                                     X14 = col_skip()), 
                    skip = ifelse(tail(readLines(filname), 1) == "", 1, 0) , 
                    trim_ws = TRUE)
  
  colnames(res1) <- c("k", "me", "mte"
                     ,"dc" , "rc" ,"xc" 
                     , "cc"  , "bc"  
                     , "ec" , "kc" ,  "pc"  
                     , "hc" , "pr")
  
  #check if all is OK
  if ( table(res1$k)[1]==2809 & sum(table(res1$me)==53)==53 & sum(table(res1$mte)==53)==53 ) 
      {
      res1 <- data.table(res1)
      res <- rbind(res, res1)
      res1 <- NULL
      }
}

#checking if everything with OK, 407 NAs at bc
summary(res)

#some bc values are NA-s, calculate them and we consider them as 0.5
sum(is.na(res$bc))
res$bc[is.na(res$bc)] <- 0.5

#05. ANALYSES -- CREATING the qetb table 95% p error tau b values
#05. ANALYSES -- CREATING the qetb table 95% p error tau b values
head(res)

#saveRDS(res, file="resraw_23thMay2022_k128")
#res<-readRDS(file="resraw_22thMay2022")


#using etb(i,k,me,mte) calculate 95% quantiles of qetbi(me,mte)
qetb <- res %>%
  group_by(me, mte) %>%
  summarise(across(everything(), ~ quantile(., probs = 0.95, na.rm = TRUE))) 

qetb <- qetb %>%
  rowwise %>%
  mutate (mer=me/ec1) %>%
  rowwise %>%
  mutate (mter=mte/ec1)

#saveRDS(qetb, file="qetb_23May2023_k128.RDS")  

#Figure1 - colouring for visualization
qetb$c_dc <- cut(qetb$dc, 
                       breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                       labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

qetb$c_rc <- cut(qetb$rc, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

qetb$c_xc <- cut(qetb$xc, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

#Figure1 - colouring for visualization
qetb$c_cc <- cut(qetb$cc, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

#Figure1 - colouring for visualization
qetb$c_bc <- cut(qetb$bc, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

#Figure1 - colouring for visualization
qetb$c_ec <- cut(qetb$ec, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

#Figure1 - colouring for visualization
qetb$c_kc <- cut(qetb$kc, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

#Figure1 - colouring for visualization
qetb$c_pc <- cut(qetb$pc, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

#Figure1 - colouring for visualization
qetb$c_hc <- cut(qetb$hc, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

#Figure1 - colouring for visualization
qetb$c_pr <- cut(qetb$pr, 
                 breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                 labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))


#check is there is na or nan
sum(is.na(qetb))
sum(is.nan(qetb$bc))

#Figure 1. Boxplot
qetb_long <-  qetb[c(1,2,4,5,6,7,8,9,10,11,12,13)] %>%
  tidyr::gather(cmeas, etb, dc:pr)

#colouring for visualization
qetb_long$c <- cut(qetb_long$etb, 
                   breaks = c(-Inf, 0.05, 0.1, 0.2, 0.3, Inf), 
                   labels = c("minimal error", "small error", "moderate error", "large error", "excessive error"))

ggplot(qetb_long[(qetb_long$me + qetb_long$mte > 10) & (qetb_long$me + qetb_long$mte < 15) ,]
       , aes(x = reorder(cmeas, etb, FUN=median), y = etb)) +
  geom_boxplot(aes(fill = cmeas)) +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -1, hjust = 0.5, size = 4, color = "black") +
  xlab("Centrality measures") +
  ylab("95% percentile of the proportion of discordant pairs (qetb)") +
  ggtitle("Sensitivity in the 3.9-5.9% noise area") +
  scale_fill_manual(values = c("yellow","orange", "lightgreen", "yellow", "yellow", 
                               "yellow", "orange","yellow", "lightgreen", "orange")) +
  guides(colors="none", fill="none")


ggplot(qetb_long[(qetb_long$me + qetb_long$mte > 22) & (qetb_long$me + qetb_long$mte < 29) ,]
       , aes(x = reorder(cmeas, etb, FUN=median), y = etb)) +
  geom_boxplot(aes(fill = cmeas)) +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -1, hjust = 0.5, size = 4, color = "black") +
  xlab("Centrality measures") +
  ylab("95% percentile of the proportion of discordant pairs (qetb)") +
  ggtitle("Sensitivity in the 8.6-11.4% noise area") +
  scale_fill_manual(values = c("yellow","orange", "lightgreen", "yellow", "yellow", 
                               "yellow", "orange","yellow", "lightgreen", "orange")) +
  guides(colors="none", fill="none")



#figure2a_degree
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_dc)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=1, Degree centrality") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(dc*100, 1))), color = "black", size = 2) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2b_Reach(k=2)
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_rc)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=2, Reach centrality (k=2)") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(rc*100, 1))), color = "black", size = 3) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2c_Eccentricity
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_ec)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=3, Eccentricity") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(ec*100, 1))), color = "black", size = 3) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2d_Closeness
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_cc)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=4, Closeness centrality") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(cc*100, 1))), color = "black", size = 2) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2e_Betwenness
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_bc)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=5, Betwenness centrality") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(bc*100, 1))), color = "black", size = 3) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2f_Eigenvector
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_ec)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=6, Eigenvector centrality") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(ec*100, 1))), color = "black", size = 3) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2g_Kleinberg_authority
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_kc)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=7, Kleinberg authority centrality (l=0.9)") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(kc*100, 1))), color = "black", size = 3) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2h_Power_centrality
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_pc)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=8, Power centrality") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(pc*100, 1))), color = "black", size = 2) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2i_Harmonic_centrality
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_hc)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=9, Harmonic centrality") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(hc*100, 1))), color = "black", size = 3) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))

#figure2j_Page rank
ggplot(qetb, aes(x = (mer*100), y = (mter*100), fill = c_pr)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly removed edges (percentage)", y ="randomly added edges (percentage) ", fill="qetb",
       title = "i=10, Page rank centrality") +
  xlim(-0.5,5.7) +
  ylim(-0.5,5.7) +
  geom_text(aes(label = sprintf("%.1f", round(pr*100, 1))), color = "black", size = 2) + 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "yellow", "orange", "red", "darkred"))


##we cc, pc and xc because they are very noisy!
res <- subset(res, select = -c(cc, pc, xc))

#06. ANALYSES -- APPROXIMATIONS
#06. ANALYSES -- APPRORXIMATIONS
#res<-readRDS(file="resraw_22thMay2022")

fexp <- function(x1, x2, b1, b2) {(1 - (exp((-b1 * x1) + (-b2*x2) ) )) }
flog <- function(x1, x2, b1, b2) {2*((1 / (1 + exp( (-b1 * x1) -  (b2 * x2) )))-0.5)}
fhypt <- function(x1, x2, b1, b2) {1 * (tanh( (b1 * x1) + (b2*x2) ) ) }
fgud <- function (x1, x2, b1, b2) { sqrt(2) * atan(tanh(b1*x1 + b2*x2)) }
fgpd <- function (x1, x2, b1, b2, c) {1 * ( (b1*x1 + b2*x2)/ (1+((b1*x1 + b2*x2)^c))^(1/c)) }
fgef <- function (x1, x2, b1, b2) 
  {ga1 = 0.278393 
  ga2 = 0.230389
  ga3 = 0.000972
  ga4 = 0.078108
  1 - (1 / ((1+ ga1*(b1*x1 + b2*x2) + ga2*(b1*x1 + b2*x2)^2 +  ga3*(b1*x1 + b2*x2)^3 + ga4*(b1*x1 + b2*x2)^4)^4)) }

#visual checking if my functions do what they should
tmpx <- seq(0, 0.2, by = 0.01)
plot(fexp(tmpx,tmpx,5,5)~tmpx, col="black", pch = 19, xlab = "x", ylab="y", ylim=c(0,1))
points(tmpx, flog(tmpx,tmpx,5,5), col="lightgreen", pch = 19)
points(tmpx, fhypt(tmpx,tmpx,5,5), col="blue", pch = 19)
points(tmpx,fgud(tmpx,tmpx,5,5), col="purple", pch = 19)
points(tmpx,fgpd(tmpx,tmpx,5,5,2), col="orange", pch = 19)
points(tmpx,fgef(tmpx,tmpx,5,5), col="yellow", pch = 19)



#finding the coefficients
#res is longed
head(res)
resl <- res[,c(2:10)] %>% 
  tidyr::gather(cmeas, etb, dc:pr)
head(resl)

#adding noise ratio variable
resl <- resl %>%
  rowwise %>%
  mutate (mer=me/ec1) %>%
  rowwise %>%
  mutate (mter=mte/ec1) %>%
  rowwise %>%
  mutate (noise=mer + mter)

head(resl)

#cutting it to 20% noise level
resl <- resl[(resl$noise)<0.2033,]


#plot for finding approximations
#Figure2
ggplot(resl[ (resl$cmeas=="pr") & ((resl$mer + resl$mer) <0.201) ,], aes(x = (mer + mter), y = etb, color=cmeas, fill=cmeas)) +
  geom_point(alpha=0.05) +
  stat_summary(fun.y = median, geom = "line", size = 1, color="black") +
  labs(x = " Noise ratio (proportion of removed edges and proportion of added edges) ",
       y = "Proportion of discordant pairs (qetb) [black line is the median]") +
  #geom_abline(intercept = 0, slope=1.625, linetype = "dashed") +
  guides(color="none", fill="none") 


#saveRDS(resl, file = "resl_23May2023_k128.RDS")

#nls approximations
bs_fexp <- summary(nls( (2*etb)  ~  fexp(mter,mer,b1,b2) 
             , start = list(b1 = 0, b2 = 0)
             , data = resl)
  )
b1_fexp <- bs_fexp$coefficients[1]
b2_fexp <- bs_fexp$coefficients[2]


bs_flog <- summary(nls( (2*etb)  ~  flog(mter,mer,b1,b2)
             , start = list(b1 = 0, b2 = 0)
             , data = resl,
             control = nls.control(maxiter = 100))
          )
b1_flog <- bs_flog$coefficients[1]
b2_flog <- bs_flog$coefficients[2]


bs_fhypt <- summary(nls( (2*etb)  ~  fhypt(mter,mer,b1,b2) 
             , start = list(b1 = 0, b2 = 0)
             , data = resl)
  )
b1_fhypt <- bs_fhypt$coefficients[1]
b2_fhypt <- bs_fhypt$coefficients[2]

bs_fgud <- summary(nls( (2*etb)  ~  fgud(mter,mer,b1,b2) 
             , start = list(b1 = 0, b2 = 0)
             , data = resl)
  )
b1_fgud <- bs_fgud$coefficients[1]
b2_fgud <- bs_fgud$coefficients[2]

#here the problem is complex, hence
#the first modell is for starting value approximation
bs_fgpd_tmp <- summary(nls( (2*etb) ~ c+ (b1*mter + b2*mer) / (1+(b1*mter + b2*mer))
                        , start = list(b1 = 0, b2 =  0, c = 1 )
                        , data = resl)
      )
  
bs_fgpd <- summary(nls( (2*etb)  ~  fgpd(mter,mer,b1,b2, 2) 
             , start = list(b1 = 0, b2 =  0)
             , data = resl
             )
          )
b1_fgpd <- bs_fgpd$coefficients[1]
b2_fgpd <- bs_fgpd$coefficients[2]



bs_fgef <- summary(nls( (2*etb)  ~  fgef(mter,mer,b1,b2) 
             , start = list(b1 = 0, b2 = 0)
             , data = resl)
  )
b1_fgef <- bs_fgef$coefficients[1]
b2_fgef <- bs_fgef$coefficients[2]


etb_res<- data.frame(i = c("Expd","Log","Hypt","Gud","Gpd", "Gef"), 
                    b1 = c(b1_fexp,b1_flog,b1_fhypt,b1_fgud,b1_fgpd, b1_fgef),
                    b2 = c(b2_fexp,b2_flog,b2_fhypt,b2_fgud,b2_fgpd, b2_fgef)
                  )
etb_res[,c(2,3)] <- round(etb_res[,c(2,3)],2)
View(etb_res)

#three parameter functions, the cs were set to:
c1_fgpd <- 2

#creating the prediction variables
resl$etb_pred_expd <- fexp(resl$mter, resl$mer, b1_fexp, b2_fexp) / 2
resl$etb_pred_log <- flog(resl$mter, resl$mer, b1_flog, b2_flog) / 2
resl$etb_pred_hypt <- fhypt(resl$mter, resl$mer, b1_fhypt, b2_fhypt) / 2
resl$etb_pred_gud <- fgud(resl$mter, resl$mer, b1_fgud, b2_fgud) / 2
resl$etb_pred_gpd <- fgpd(resl$mter, resl$mer, b1_fgpd, b2_fgpd, c1_fgpd) / 2
resl$etb_pred_gef <- fgef(resl$mter, resl$mer, b1_fgef, b2_fgef) / 2

#Model selection
mse_res <- NULL
#full dataset
mse_expd <- sqrt(abs(resl$etb - resl$etb_pred_expd))
mse_log <- sqrt(abs(resl$etb - resl$etb_pred_log))
mse_hypt <- sqrt(abs(resl$etb - resl$etb_pred_hypt))
mse_gud <- sqrt(abs(resl$etb - resl$etb_pred_gud))
mse_gpd <- sqrt(abs(resl$etb - resl$etb_pred_gpd))
mse_gef <- sqrt(abs(resl$etb - resl$etb_pred_gef))

#adding to a dataframe
mse_res <- data.frame(
    expd = mse_expd
    ,log = mse_log
    ,hypt = mse_hypt
    ,gud = mse_gud
    ,gpd = mse_gpd
    ,gef = mse_gef
    )

#finding the coefficients
#res is longed
head(mse_res)
mse_res_long <- mse_res %>% 
  tidyr::gather(model, msre, expd:gef)
head(mse_res_long)

ggplot(mse_res_long , aes(x = reorder(model, msre), y = msre)) +
  geom_boxplot(aes(fill = model)) +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 3)), 
               vjust = -1, hjust = 0.5, size = 3, color = "black") +
  xlab("Model") +
  ylab("Square root of absolute errors")  +
  scale_fill_manual(values = c("lightgreen","yellow", "yellow", "yellow", "yellow", "yellow")) +
  guides(colors="none", fill="none")

anovamod <- aov(msre ~ model, data = mse_res_long)
summary(anovamod)
posthoc <- TukeyHSD(anovamod, conf.level = 0.99)
filtered_posthoc <- as.data.frame(as.matrix(posthoc$`model`)) %>%
    filter(grepl("expd$", rownames(as.data.frame(as.matrix(posthoc$`model`)))))
min(filtered_posthoc$lwr)
max(filtered_posthoc$lwr)

#we went on with expd
#the reason log is worst than expd is because no +- parameters were allow for none of the models!


#going back to both dimensions just for visualization
#Figure 2a. Exp
ggplot(resl, aes(x = mer+mter, y = etb)) +
  geom_point(alpha=0.1) +
  labs(x = " Noise (proportion of removed edges + proportion of added edges)",
       y = "Raw etb values") +
  geom_point(aes (y = etb_pred_expd , color="red") ) +
  guides(color="none")

loess <- loess(etb ~ mer + mter, data = resl)

ggplot(resl, aes(x = mer+mter, y = etb)) +
  geom_point(alpha=0.1) +
  labs(x = " Proportion of removed edges + proportion of added edges",
       y = "Raw etb values") +
  geom_point(aes (y = etb_pred_log , color="red") ) +
  geom_smooth(method = "loess", se=FALSE) +
  guides(color="none")

ggplot(resl, aes(x = mer+mter, y = etb)) +
  geom_point(alpha=0.1) +
  labs(x = " Proportion of removed edges + proportion of added edges",
       y = "Raw etb values") +
  geom_point(aes (y= etb_pred_hypt , color="red") ) +
  guides(color="none")

ggplot(resl, aes(x = mer+mter, y = etb)) +
  geom_point(alpha=0.1) +
  labs(x = " Proportion of removed edges + proportion of added edges",
       y = "Raw etb values") +
  geom_point(aes (y= etb_pred_gud , color="red") ) +
  guides(color="none")

ggplot(resl[(resl$mer<0.1 & resl$ier<0.1),], aes(x = mer+ier, y = etb)) +
  geom_point(alpha=0.1) +
  labs(x = " Proportion of removed edges + proportion of added edges",
       y = "Raw etb values") +
  geom_point(aes (y= etb_pred_gpd  , color="red") ) +
  guides(color="none")

ggplot(resl, aes(x = mer+mter, y = etb)) +
  geom_point(alpha=0.1) +
  labs(x = " Proportion of removed edges + proportion of added edges",
       y = "Raw etb values") +
  geom_point(aes (y= etb_pred_gef , color="red") ) +
  guides(color="none")
 
#let's take a loot by centrality measures using exponential decrease

#adding "noise ratio" variable
res <- res %>%
  rowwise %>%
  mutate (mer=me/ec1) %>%
  rowwise %>%
  mutate (mter=mte/ec1) %>%
  rowwise %>%
  mutate (noise=mer + mter)

#only analzing 0-20% noise
res <- res[(res$noise)<0.2033,]

dcm <- summary(nls(2*dc ~ fexp(mter, mer,b1,b2), start = list(b1 = 0, b2= 0), data = res))
rcm <- summary(nls(2*rc ~ fexp(mter, mer,b1,b2), start = list(b1 = 0, b2= 0), data = res))

kcm <- summary(nls(2*kc ~ fexp(mter, mer,b1,b2), start = list(b1 = 0, b2= 0), data = res))
ecm <- summary(nls(2*ec ~ fexp(mter, mer,b1,b2), start = list(b1 = 0, b2= 0), data = res))

prm <- summary(nls(2*pr ~ fexp(mter, mer,b1,b2), start = list(b1 = 0, b2= 0), data = res))
hcm <- summary(nls(2*hc ~ fexp(mter, mer,b1,b2), start = list(b1 = 0, b2= 0), data = res))

bcm <- summary(nls(2*bc ~ fexp(mter, mer,b1,b2), start = list(b1 = 0, b2= 0), data = res))


dcb1 <- dcm$coefficients[1]
dcb2 <- dcm$coefficients[2]

rcb1 <- rcm$coefficients[1]
rcb2 <- rcm$coefficients[2]

kcb1 <- kcm$coefficients[1]
kcb2 <- kcm$coefficients[2]

ecb1 <- ecm$coefficients[1]
ecb2 <- ecm$coefficients[2]

prb1 <- prm$coefficients[1]
prb2 <- prm$coefficients[2]

hcb1 <- hcm$coefficients[1]
hcb2 <- hcm$coefficients[2]

bcb1 <- bcm$coefficients[1]
bcb2 <- bcm$coefficients[2]


nl_res<- data.frame(i = c("Degree","Reach(k=2)","Kleinberg","Eigenvector","Page rank","Harmonic","Betweenness"), 
                    b1 = c(dcb1,rcb1,kcb1,ecb1,prb1,hcb1,bcb1),
                    b2 = c(dcb2,rcb2,kcb2,ecb2,prb2,hcb2,bcb2),
                    b1xb2 = c(dcb1,rcb1,kcb1,ecb1,prb1,hcb1,bcb1) * c(dcb2,rcb2,kcb2,ecb2,prb2,hcb2,bcb2)
                  )

nl_res[,c(2,3,4)] <- round(nl_res[,c(2,3,4)],2)
View(nl_res)

res$dc_pred <-  fexp(res$mter, res$mer, dcb1, dcb2) / 2
res$rc_pred <-  fexp(res$mter, res$mer, rcb1, rcb2) / 2
res$kc_pred <-  fexp(res$mter, res$mer, kcb1, kcb2) / 2
res$ec_pred <-  fexp(res$mter, res$mer, ecb1, ecb2) / 2
res$pr_pred <-  fexp(res$mter, res$mer, prb1, prb2) / 2
res$hc_pred <-  fexp(res$mter, res$mer, hcb1, hcb2) / 2
res$bc_pred <-  fexp(res$mter, res$mer, bcb1, bcb2) / 2


#full dataset
mse_dc <- sqrt(abs(res$dc - res$dc_pred))
mse_rc <- sqrt(abs(res$rc - res$rc_pred))
mse_kc <- sqrt(abs(res$kc - res$kc_pred))
mse_ec <- sqrt(abs(res$ec - res$ec_pred))
mse_pr <- sqrt(abs(res$pr - res$pr_pred))
mse_hc <- sqrt(abs(res$hc - res$hc_pred))
mse_bc <- sqrt(abs(res$bc - res$bc_pred))

summary(mse_dc)
summary(mse_rc)
summary(mse_kc)
summary(mse_ec)
summary(mse_hc)
summary(mse_bc)
summary(mse_pr)
