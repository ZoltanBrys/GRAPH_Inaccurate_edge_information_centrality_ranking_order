#Author1=ORCID:0000-0002-3324-2255
#Author2=ORCID:0000-0002-1319-0785
#Title=Effect of inaccurate edge information on centrality ranking error in classroom networks
#ver=3.221
#date=10th of May, 2025

#00.ENVIRONMENT and PACKAGES
#00.ENVIRONMENT and PACKAGES
rm(list = ls())
library("igraph")
library("readr")
library("dplyr")
library("tidyr")
library("googledrive")
library("ggplot2")


#calculating confidence intervall for kendall-tau-b
#fieller correction:
#Fieller, E. C., Hartley, H. O., & Pearson, E. S. (1957). 
#Tests for rank correlation coefficients: I. 
#Biometrika, 44, 470–481.
ktau.ci <- function(tau, N, conf.level = 0.95, correct='fieller') {
  if(correct=='none') tau.se <- 1/(N - 3)^0.5
  if(correct=='fieller') tau.se <- (0.437/(N - 4))^0.5
  moe <- qnorm(1 - (1 - conf.level)/2) * tau.se
  zu <- atanh(tau) + moe
  zl <- atanh(tau) - moe
  tanh(c(zl, zu))
}

#01.FOR PRESENTATION and DOUBLE-CHECKING HOW R CALCULATES KENDALL-TAU-B
#01.FOR PRESENTATION and DOUBLE-CHECKING HOW R CALCULATES KENDALL-TAU-B
refc<-c(1:141)
ref1a<-c(c(1:118),c(141:119)) 
cor(refc, ref1a , use = "complete.obs", method = "kendall")
cor(refc, ref1a , use = "complete.obs", method = "kendall")*(141*140)
#good, it calculates the way it shall

#02.RAW DATA and PREPATATION of the DATA
#02.RAW DATA and PREPATATION of the DATA

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


#REATING the "real" graph, rg
#rg = real graph
rg <- graph_from_data_frame(cnss, directed = FALSE)
rg <- simplify(rg, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = igraph_opt("edge.attr.comb"))

#taking the biggest component from centrality measures it is more comprehensible if the is one component
comp <- components(rg)
rg <- decompose(rg, max.comps = 1, min.vertices = max(comp$csize))[[1]] #selecting the biggest component
rm(comp)

#for the presentation only
#imap1 <- cluster_infomap(rg)
#plot(imap1 , rg, edge.color = "transparent")

#basic parameters of the real graph
vcount(rg)
ecount(rg)
plot(rg, layout=layout_with_fr(rg), vertex.size=5, vertex.label.cex=0.3)
table(degree(rg))
hist(degree(rg))
median(degree(rg))
mean(degree(rg))

#degree sequence is to be saved for simulating alike graphs for calculating baseline error
rg_degree_seq <- degree(rg)

#ground truth centrality parameters for each node
rg_cc <- closeness(rg)
rg_bc <- betweenness(rg)
rg_ec <- eigen_centrality(rg)$vector
rg_hc <- harmonic_centrality(rg)
rg_pr <- page_rank(rg)$vector


rg_ec1 <- ecount(rg) #edge count of the real graph
rg_nc1 <- vcount(rg) #node count of the real graph

rg_ec1

#03.SIMULATION of the LIFE-LIKE ZONE DISTORTED GRAPHS and tau-b calculation between original and distorted graphs
#03.SIMULATION of the LIFE-LIKE ZONE DISTORTED GRAPHS and tau-b calculation between original and distorted graphs
res<-NULL #results table

for (k in 1:16)
 for (me in 0 : rg_ec1 )
  for (ie in 0 : rg_ec1 )
    {
    
    if (me==0) {sg2<-rg} else #if edges are removed
        {
         deids <- sample(x=c(1:rg_ec1), size=me, replace=FALSE) #ids of edges to delete
         sg1 <- delete_edges(rg, deids) #delete
         sg2 <- sg1 #pass on
        }
    
    #if the number of inaccurate edges are more then all of the edges, than max it.
    if (ecount(sg2)>ie) { ie2 <- ie } else { ie2 <- (ecount(sg2)) } 
      
    if (!ie2==0) {
                edges_to_modify <- sample(E(sg2), size = ie2) 
                sg2 <- delete_edges(sg2, edges_to_modify) #edges deleted from the graph
                
                  for (c in 1:length(edges_to_modify) ) #now let' add new ones
                      {
                       nodes <- sample(V(sg2), size = 2)
                       while ( are.connected(sg2, nodes[1], nodes[2]) ) {nodes <- sample(V(sg2), size = 2) } 
                       sg2 <- add_edges(sg2, c(nodes[1], nodes[2]) )
                      }
                  }
    
    #now do the analysis
    #sg is the sampled, moisy, distorted graph
    sg_cc <- try(closeness(sg2))
    sg_bc <- try(betweenness(sg2))
    sg_ec <- try(eigen_centrality(sg2)$vector)
    sg_hc <- try(harmonic_centrality(sg2))
    sg_pr <- try(page_rank(sg2)$vector)
    
    
    #calculating the kendall tau based error the two full graphs
    
    tmp_cc <- 0
    tmp_bc <- 0
    tmp_ec <- 0
    tmp_hc <- 0
    tmp_pr <- 0
    
    if (sum(is.nan(sg_cc))<140) 
        {
      tmp_cc <- cor(rg_cc, sg_cc,use = "complete.obs", method = "kendall")
      tmp_bc <- cor(rg_bc, sg_bc,use = "complete.obs", method = "kendall")
      tmp_ec <- cor(rg_ec, sg_ec,use = "complete.obs", method = "kendall")
      tmp_hc <- cor(rg_hc, sg_hc,use = "complete.obs", method = "kendall")
      tmp_pr <- cor(rg_pr, sg_pr,use = "complete.obs", method = "kendall")
        }
    
    cc_er_ke <- 1 - tmp_cc
    bc_er_ke <- 1 - tmp_bc
    ec_er_ke <- 1 - tmp_ec
    hc_er_ke <- 1 - tmp_hc
    pr_er_ke <- 1 - tmp_pr
    
    cc_er_ke_l <- 1 - ktau.ci(tmp_cc, 141)[1]
    bc_er_ke_l <- 1 - ktau.ci(tmp_bc, 141)[1]
    ec_er_ke_l <- 1 - ktau.ci(tmp_ec, 141)[1]
    hc_er_ke_l <- 1 - ktau.ci(tmp_hc, 141)[1]
    pr_er_ke_l <- 1 - ktau.ci(tmp_pr, 141)[1]
    
    if (is.nan(cc_er_ke_l)) {cc_er_ke_l<-1}
    if (is.nan(cc_er_ke_l)) {cc_er_ke_l<-1}
    if (is.nan(cc_er_ke_l)) {cc_er_ke_l<-1}
    if (is.nan(cc_er_ke_l)) {cc_er_ke_l<-1}
    if (is.nan(cc_er_ke_l)) {cc_er_ke_l<-1}
    
      res1 <- data.frame(
                     k=k
                    ,me=me
                    ,ie=ie
                    ,mer =me/rg_ec1 # missing edges ratio
                    ,ier =ie/rg_ec1 # inaccurate edges ratio
                    ,cc_er_ke
                    ,bc_er_ke
                    ,ec_er_ke
                    ,hc_er_ke
                    ,pr_er_ke
                    ,cc_er_ke_l
                    ,bc_er_ke_l
                    ,ec_er_ke_l
                    ,hc_er_ke_l
                    ,pr_er_ke_l
                    ,tim=as.character(Sys.time())
                    )
    
    res<-rbind(res,res1)
    
    cat("k=", k, " missing edges=", me, " inaccurate edges=", ie, "\n")   
    
    Sys.sleep(0.02)

} 

saveRDS(res, file="EUSN_20230512_k16.RDS") #k=256, standard error zone, te included.

#CHECKING IF THE SIMULATION WAS OK
table(res$k) #it shall be 1-100 and constant 16384
table(res$me) #it shall be constant 128
table(res$ie) #it shall be constant 128

#05.ANALYYING AND VISUALIZING the RESULTS of RAW tau-beta and baselines
#05.ANALYYING AND VISUALIZING the RESULTS of RAW tau-beta and baselines
#library(psych) #for anova
#library(car) #for anova
#library(effectsize) #for effect size

#first one-way-anova for all rbes by centrality measure, ANOVA
res_long <- (((res[,c(2,3,6,7,8,9,10)]))) %>% 
  tidyr::gather(cmeas, taub, cc_er_ke:pr_er_ke)

oneway.test(taub~cmeas, data = res_long, var.equal = FALSE)
kruskal.test(taub~cmeas, data = res_long)
oneway_lm <- lm(taub~cmeas, data = res_long)
summary(oneway_lm)

512*0.05

res <- res %>%
  rowwise() %>%
  mutate(min_tbe = min(cc_er_ke, bc_er_ke, ec_er_ke, hc_er_ke, pr_er_ke))

ggplot(res, aes(x = me, y = ie, fill = min_tbe)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2(low = "white", mid = "red", high = "darkred",
                       midpoint = 0.5, limits = c(0, 2)) +
  labs(x = "randomly  removed edges", y = "randomly added edges ", fill="τbe", 
       title = "Minimal ranking errors of the five centrality measure") + 
  theme(plot.title = element_text(size=10))

#only the first 8x8
res$cg_min_tbe <- cut(res$min_tbe, breaks = c(Inf, 0.1, 0.3, -Inf), 
                      labels = c("small error","moderate error", "large error"))
ggplot(res, aes(x = me, y = ie, fill = cg_min_tbe)) +
  geom_tile() +
  coord_equal()  +
  labs(x = "randomly  removed edges (0-256)", y ="randomly added edges (0-256) ", fill="τbe",
       title = "Minimal ranking errors of the five centrality measure") +
  xlim(-0.5,17.5) +
  ylim(-0.5,17.5) +
  geom_text(aes(label = sprintf("%.1f", round(min_tbe, 1))), color = "black", size = 3) + 
#  geom_text(aes(label = round(min_tbe, 1)), color = "black", size = 3)+ 
  theme(plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("green", "orange", "darkred"))

#comparing performance in the zone
ggplot(res_long[ (res_long$ie + res_long$me) < 17 ,], aes(x = cmeas, y = taub)) +
  geom_boxplot(outlier.shape = 19, notch = FALSE) +
  stat_summary(fun.y = median, geom = "point", shape = 18, size = 3, color = "red") +
  stat_summary(fun.y = median, geom = "text", vjust=-1, color="red",
               aes(label=sprintf("%.2f",..y..)),
               size = 3) +
  scale_x_discrete(labels = c(cc_er_ke = "Closeness", 
                              bc_er_ke = "Betweenness",
                              ec_er_ke = "Eigenvector",
                              hc_er_ke = "Harmonic",
                              pr_er_ke = "Page rank"))

#Defining noise ratio and adding to res
res <- res %>%
  rowwise() %>%
  mutate(noise = (ie + me) ) %>%
  rowwise() %>%
  mutate(noiser = (noise/512) )

summary(res$noise)
summary(res$noiser)

res_long <- (((res[res$noise<257,c(19,2,3,6,7,8,9,10)]))) %>% 
  tidyr::gather(cmeas, taub, cc_er_ke:pr_er_ke)

ggplot(res_long, aes(x = round((noiser),2)*100, group=cmeas, color=cmeas, y = taub)) +
  geom_line(alpha=0) +
  stat_summary(fun.y = median, geom = "point", shape = 18, size = 3) +
  labs(x = "Noise ratio (0-100%)", y = "tbe") +
  geom_hline(yintercept = 0.1, color="orange") +
  geom_hline(yintercept = 0.3, color="red") +
  scale_color_discrete(labels = c(cc_er_ke = "Closeness", 
                              bc_er_ke = "Betweenness",
                              ec_er_ke = "Eigenvector",
                              hc_er_ke = "Harmonic",
                              pr_er_ke = "Page rank")) +
  coord_cartesian(ylim=c(0,1.1))

#beta regression
library(betareg)

res$n_cc_er_ke <- (res$cc_er_ke/2)
res$n_bc_er_ke <- (res$bc_er_ke/2)
res$n_ec_er_ke <- (res$ec_er_ke/2)
res$n_hc_er_ke <- (res$hc_er_ke/2)
res$n_pr_er_ke <- (res$pr_er_ke/2)

hist(res$n_bc_er_ke)
summary(res$n_cc_er_ke)

betareg(n_cc_er_ke ~ mer + ier , data=res[res$noise<512,])
betareg(n_bc_er_ke ~ mer + ier , data=res[res$noise<128,])
betareg(n_ec_er_ke ~ mer + ier , data=res)
betareg(n_hc_er_ke ~ mer + ier , data=res)
betareg(n_pr_er_ke ~ mer + ier , data=res)

#For the EUSN presentation
#For the EUSN presentation
refc<-c(1:10)
ref1a<-c(1,2,3,4,5,6,7,8,10,9)
ref1b<-c(10,2,3,4,5,6,7,8,9,1)
ref1c<-c(5,6,4,7,3,8,9,2,1,10)
ref1e<-c(9,10,8,7,6,5,4,3,2,1)
ref1g<-c(10,9,8,7,6,5,4,3,2,1)



cor(refc, ref1a , use = "complete.obs", method = "kendall")*90
cor(refc, ref1b , use = "complete.obs", method = "kendall")*90
cor(refc, ref1c , use = "complete.obs", method = "kendall")*90
cor(refc, ref1e , use = "complete.obs", method = "kendall")
cor(refc, ref1g , use = "complete.obs", method = "kendall")


figure1<-data.frame(
  n1=c(       "Sandra",
           "John", 
           "Michael", 
           "Mary",
           "Sandra", 
           "Michael", 
            "James",
            "James",
           "Sandra",
           "Sandra"),
   n2=c( 
               "John", 
               "Patrick", 
               "Mary", 
               "Patrick",
               "Karen", 
               "John", 
               "John",
               "Patrick",
               "Jennifer",
               "David"))
  
#figure 1 graph
f1g <- graph_from_data_frame(figure1, directed = FALSE)
hcf1g <- rank(harmonic_centrality(f1g))

gl <- layout_with_fr(f1g)

par(mar = c(0, 0, 0, 0))
plot(f1g, layout=gl, vertex.size=20, vertex.label.cex=1.2, vertex.label = paste0(V(f1g)$name))

par(mar = c(0, 0, 0, 0))
plot(f1g, layout=gl, vertex.size=20, vertex.label.cex=1.2, 
     vertex.label = paste0(V(f1g)$name, "   (", hcf1g, ")") )

f1g<-delete_edges(f1g,2)
nhcf1g <- rank(harmonic_centrality(f1g))
par(mar = c(0, 0, 0, 0))

par(mar = c(0, 0, 0, 0))
plot(f1g, layout=gl, vertex.size=20, vertex.label.cex=1.2, vertex.label = paste0(V(f1g)$name))

plot(f1g, layout=gl, vertex.size=20, vertex.label.cex=1.2, 
     vertex.label = paste0(V(f1g)$name, " (", nhcf1g, ")") )

dif<-data.frame(name=as.character(V(f1g)$name), realhc=hcf1g, disthc=nhcf1g )
dif$dif<-dif$realhc-dif$disthc
dif
cor.test()
