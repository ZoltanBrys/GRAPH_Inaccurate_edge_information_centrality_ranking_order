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

#overlap function is to estimate the overlap between nodes
overlap1 <- function(x, y) {
  if (length(x)==length(x))
    {
  common_names <- intersect(names(x), names(y))
  
  return( length(common_names) / length(x) ) 
  }
}

#01.FOR PRESENTATION and DOUBLE-CHECKING HOW R CALCULATED KENDALL-TAU-B
#01.FOR PRESENTATION and DOUBLE-CHECKING HOW R CALCULATED KENDALL-TAU-B
refc<-c(1:5)
ref1a<-c(1,2,3,4,5)
ref1b<-c(1,2,3,5,4)
ref1c<-c(5,2,3,5,1)
ref1d<-c(2,1,3,5,4)
ref1e<-c(2,5,1,3,4)

cor(refc, ref1a , use = "complete.obs", method = "kendall")
cor(refc, ref1b , use = "complete.obs", method = "kendall")
cor(refc, ref1c , use = "complete.obs", method = "kendall")
cor(refc, ref1d , use = "complete.obs", method = "kendall")
cor(refc, ref1e , use = "complete.obs", method = "kendall")

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


#01. CREATING the "real" graph, rg
#01. CREATING the "real" graph, rg
#rg = real graph
rg <- graph_from_data_frame(cnss, directed = FALSE)
rg <- simplify(rg, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = igraph_opt("edge.attr.comb"))

#taking the biggest component from centrality measures it is more comprehensible if the is one component
comp <- components(rg)
rg <- decompose(rg, max.comps = 1, min.vertices = max(comp$csize))[[1]] #selecting the biggest component
rm(comp)

#basic parameters of the real graph
vcount(rg)
ecount(rg)
plot(rg, layout=layout_with_fr(rg))
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


#mean(degree(rg))
#mean degree is 3.3
#on an average day 8.2% of pupils missing 5%-10%
#on an average day from this class
# 185 * 8.2% missing from 185 * 5% - 185 *10%

#MAIN CYCLE
#k = repeat k times, because it is stochastic
#me = missing edge
#ie = inaccurate edge
#te = top edges 1-100%

rg_ec1 <- ecount(rg) #edge count of the real graph
rg_nc1 <- vcount(rg) #node count of the real graph

#real-life zone calculations
round(rg_ec1*0.04) 
round(rg_ec1*0.1)

round(rg_ec1*0.01)
round(rg_ec1*0.03)


#03.SIMULATION of the LIFE-LIKE ZONE DISTORTED GRAPHS and tau-b calculation between original and distorted graphs
#03.SIMULATION of the LIFE-LIKE ZONE DISTORTED GRAPHS and tau-b calculation between original and distorted graphs
res<-NULL #results table

for (k in 1:100)
 for (me in 10 : 26 )
  for (ie in 3 : 8 )
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
    #sg is sampled, distorted graph
    sg_cc <- closeness(sg2)
    sg_bc <- betweenness(sg2)
    sg_ec <- eigen_centrality(sg2)$vector
    sg_hc <- harmonic_centrality(sg2)
    sg_pr <- page_rank(sg2)$vector
    
    
    nesg2 <- vcount(sg2) #number of nodes sg2
    
            for (te in 14:nesg2) #calculating tau for top10%-100% of edges, 14 is set to achieve comparability of semi-def tau
            {
    
    #CALCULATING raw tau-b ERROR MEAUSRES, ordered and first selected n
    cc_er_ke <- cor(rg_cc[order(rg_cc, decreasing=TRUE)][1:te], sg_cc[order(rg_cc, decreasing=TRUE)][1:te]
                    , use = "complete.obs", method = "kendall")
    bc_er_ke <- cor(rg_bc[order(rg_bc, decreasing=TRUE)][1:te], sg_bc[order(rg_bc, decreasing=TRUE)][1:te]
                    , use = "complete.obs", method = "kendall")
    ec_er_ke <- cor(rg_ec[order(rg_ec, decreasing=TRUE)][1:te], sg_ec[order(rg_ec, decreasing=TRUE)][1:te]
                    , use = "complete.obs", method = "kendall")
    hc_er_ke <- cor(rg_hc[order(rg_hc, decreasing=TRUE)][1:te], sg_hc[order(rg_hc, decreasing=TRUE)][1:te]
                    , use = "complete.obs", method = "kendall")
    pr_er_ke <- cor(rg_pr[order(rg_pr, decreasing=TRUE)][1:te], sg_pr[order(rg_pr, decreasing=TRUE)][1:te]
                    , use = "complete.obs", method = "kendall")
    
    #CALCULATING overlap for comparison
    cc_er_ol <- overlap1(rg_cc[order(rg_cc, decreasing=TRUE)][1:te], sg_cc[order(rg_cc, decreasing=TRUE)][1:te])
    bc_er_ol <- overlap1(rg_bc[order(rg_bc, decreasing=TRUE)][1:te], sg_bc[order(rg_bc, decreasing=TRUE)][1:te])
    ec_er_ol <- overlap1(rg_ec[order(rg_ec, decreasing=TRUE)][1:te], sg_ec[order(rg_ec, decreasing=TRUE)][1:te])
    hc_er_ol <- overlap1(rg_hc[order(rg_hc, decreasing=TRUE)][1:te], sg_hc[order(rg_hc, decreasing=TRUE)][1:te])
    pr_er_ol <- overlap1(rg_pr[order(rg_pr, decreasing=TRUE)][1:te], sg_pr[order(rg_pr, decreasing=TRUE)][1:te])
    
    res1 <- data.frame(
                     k=k
                    ,me=me
                    ,ie=ie
                    ,mer =me/rg_ec1 # missing edges ratio
                    ,ier =ie/rg_ec1 # inaccurate edges ratio
                    ,te =te
                    ,ter = te / nesg2
                    ,nedt = nesg2 #number of nodes of the distorted graph
                    ,cc_er_ke
                    ,bc_er_ke
                    ,ec_er_ke
                    ,hc_er_ke
                    ,pr_er_ke
                    ,cc_er_ol
                    ,bc_er_ol
                    ,ec_er_ol
                    ,hc_er_ol
                    ,pr_er_ol
                    ,tim=as.character(Sys.time())
                    )
    
    res<-rbind(res,res1)
    
            }
    
    cat("k=", k, " missing edges=", me, " inaccurate edges=", ie, "\n")   
    
    Sys.sleep(0.01)
    } 

saveRDS(res, file="EUSN_2023_k100_withoverlap.RDS") #k=10, standard error zone, te included.

#04A. CALCULATING BASELINE ERRORS - Erdos Renyi
#04A. CALCULATING BASELINE ERRORS - Erdos Renyi
be_sdt_res <- NULL #baseline error for each ter and each centrality measure
be_sd_res <- NULL 

for (k in 1:100)
{
  
  rag <- erdos.renyi.game(n=vcount(rg), p.or.m=ecount(rg), type = "gnm") #random Erdos Renyi graph
  #rag1 <- sample_pa(n=vcount(rg), m=2, p=1, directed=FALSE) #creating a random graph with preferential attachment
  #shuffled_rg_degree_seq <- sample(rg_degree_seq, length(rg_degree_seq), replace = FALSE)
  #rag3 <- sample_degseq(out.deg = shuffled_rg_degree_seq, method = "simple.no.multiple.uniform")
  #rag4 <- sample_degseq(out.deg = rg_degree_seq, method = "simple.no.multiple.uniform")
  
  
  #making the number of edges equal
  if (ecount(rag)>ecount(rg))
  {
    rag <- delete.edges(rag, sample(E(rag), ecount(rag)-ecount(rg)))
  }
  
  #saving the parameters of a random graph
  sg_cc <- closeness(rag) 
  sg_bc <- betweenness(rag)
  sg_ec <- eigen_centrality(rag)$vector
  sg_hc <- harmonic_centrality(rag)
  sg_pr <- page_rank(rag)$vector
  
  nesg2 <- vcount(rag) #number of nodes random PA graph, it has to be 141
  
  for (te in 14:nesg2) #calculating tau for top1-100% of edges!
  {
    
    #CALCULATING ERROR MEAUSRES, ordered and first selected n
    cc_er_ke <- cor(rg_cc[order(rg_cc, decreasing=TRUE)][1:te], sg_cc[order(rg_cc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    bc_er_ke <- cor(rg_bc[order(rg_bc, decreasing=TRUE)][1:te], sg_bc[order(rg_bc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    ec_er_ke <- cor(rg_ec[order(rg_ec, decreasing=TRUE)][1:te], sg_ec[order(rg_ec, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    hc_er_ke <- cor(rg_hc[order(rg_hc, decreasing=TRUE)][1:te], sg_hc[order(rg_hc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    pr_er_ke <- cor(rg_pr[order(rg_pr, decreasing=TRUE)][1:te], sg_pr[order(rg_pr, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    
    res1 <- data.frame(
      k=k
      ,te =te
      ,ter = te / nesg2
      ,nedt = nesg2 #number of nodes of the random graph, has to be 141
      ,cc_er_ke
      ,bc_er_ke
      ,ec_er_ke
      ,hc_er_ke
      ,pr_er_ke
      ,tim=as.character(Sys.time())
    )
    
    be_sd_res<-rbind(be_sd_res,res1)
    
  }
  
  cat("k=", k, "\n")   
  
  Sys.sleep(0.01)
} 

colnames(be_sd_res)

#results important to note, that it start with node=3
be_sdt_res <-be_sd_res %>%
  group_by(ter) %>%
  summarize(med_cc_er_ke = median(cc_er_ke),
            med_bc_er_ke = median(bc_er_ke),
            med_ec_er_ke = median(ec_er_ke),
            med_hc_er_ke = median(hc_er_ke),
            med_pr_er_ke = median(pr_er_ke))

boxplot(be_sdt_res[be_sdt_res$ter==1, 2:6])
boxplot(be_sdt_res[be_sdt_res$ter<1, 2:6])
summary(be_sdt_res[,2:6], na.rm=TRUE)
baseline_error_erdos_renyi <- 0 # of course +-0.01

#04B. CALCULATING BASELINE ERRORS - Barabasi-Albert 
#04B. CALCULATING BASELINE ERRORS - Barabasi-Albert 
be_sdt_res <- NULL #baseline error for each ter and each centrality measure
be_sd_res <- NULL 

for (k in 1:100)
{
  
  #rag1 <- erdos.renyi.game(n=vcount(rg), p.or.m=ecount(rg), method=("gnm")) #random Erdos Renyi graph
  rag <- sample_pa(n=vcount(rg), m=2, p=1, directed=FALSE) #creating a random graph with preferential attachment
  #shuffled_rg_degree_seq <- sample(rg_degree_seq, length(rg_degree_seq), replace = FALSE)
  #rag3 <- sample_degseq(out.deg = shuffled_rg_degree_seq, method = "simple.no.multiple.uniform")
  #rag4 <- sample_degseq(out.deg = rg_degree_seq, method = "simple.no.multiple.uniform")
  
  
  #making the number of edges equal
  if (ecount(rag)>ecount(rg))
  {
    rag <- delete.edges(rag, sample(E(rag), ecount(rag)-ecount(rg)))
  }
  
  #saving the parameters of a random graph
  sg_cc <- closeness(rag) 
  sg_bc <- betweenness(rag)
  sg_ec <- eigen_centrality(rag)$vector
  sg_hc <- harmonic_centrality(rag)
  sg_pr <- page_rank(rag)$vector
  
  nesg2 <- vcount(rag) #number of nodes random PA graph, it has to be 141
  
  for (te in 14:nesg2) #calculating tau for top1-100% of edges!
  {
    
    #CALCULATING ERROR MEAUSRES, ordered and first selected n
    cc_er_ke <- cor(rg_cc[order(rg_cc, decreasing=TRUE)][1:te], sg_cc[order(rg_cc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    bc_er_ke <- cor(rg_bc[order(rg_bc, decreasing=TRUE)][1:te], sg_bc[order(rg_bc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    ec_er_ke <- cor(rg_ec[order(rg_ec, decreasing=TRUE)][1:te], sg_ec[order(rg_ec, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    hc_er_ke <- cor(rg_hc[order(rg_hc, decreasing=TRUE)][1:te], sg_hc[order(rg_hc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    pr_er_ke <- cor(rg_pr[order(rg_pr, decreasing=TRUE)][1:te], sg_pr[order(rg_pr, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    
    res1 <- data.frame(
      k=k
      ,te =te
      ,ter = te / nesg2
      ,nedt = nesg2 #number of nodes of the random graph, has to be 141
      ,cc_er_ke
      ,bc_er_ke
      ,ec_er_ke
      ,hc_er_ke
      ,pr_er_ke
      ,tim=as.character(Sys.time())
    )
    
    be_sd_res<-rbind(be_sd_res,res1)
    
  }
  
  cat("k=", k, "\n")   
  
  Sys.sleep(0.01)
} 

colnames(be_sd_res)

#results important to note, that it start with node=14
be_sdt_res <-be_sd_res %>%
  group_by(ter) %>%
  summarize(med_cc_er_ke = median(cc_er_ke),
            med_bc_er_ke = median(bc_er_ke),
            med_ec_er_ke = median(ec_er_ke),
            med_hc_er_ke = median(hc_er_ke),
            med_pr_er_ke = median(pr_er_ke))

boxplot(be_sdt_res[be_sdt_res$ter==1, 2:6])
boxplot(be_sdt_res[be_sdt_res$ter<1, 2:6])
summary(be_sdt_res[,2:6], na.rm=TRUE)

baseline_error_barabasi <- 0 # of course +-0.05

#04C. CALCULATING BASELINE ERRORS - Degree sequence kept in a random order
#04C. CALCULATING BASELINE ERRORS - Degree sequence kept in a random order
be_sdt_res <- NULL #baseline error for each ter and each centrality measure
be_sd_res <- NULL 

for (k in 1:100)
{
  
  #rag1 <- erdos.renyi.game(n=vcount(rg), p.or.m=ecount(rg), method=("gnm")) #random Erdos Renyi graph
  #rag <- sample_pa(n=vcount(rg), m=2, p=1, directed=FALSE) #creating a random graph with preferential attachment
  shuffled_rg_degree_seq <- sample(rg_degree_seq, length(rg_degree_seq), replace = FALSE)
  rag <- sample_degseq(out.deg = shuffled_rg_degree_seq, method = "simple.no.multiple.uniform")
  #rag4 <- sample_degseq(out.deg = rg_degree_seq, method = "simple.no.multiple.uniform")
  
  
  #making the number of edges equal
  if (ecount(rag)>ecount(rg))
  {
    rag <- delete.edges(rag, sample(E(rag), ecount(rag)-ecount(rg)))
  }
  
  #saving the parameters of a random graph
  sg_cc <- closeness(rag) 
  sg_bc <- betweenness(rag)
  sg_ec <- eigen_centrality(rag)$vector
  sg_hc <- harmonic_centrality(rag)
  sg_pr <- page_rank(rag)$vector
  
  nesg2 <- vcount(rag) #number of nodes random PA graph, it has to be 141
  
  for (te in 14:nesg2) #calculating tau for top1-100% of edges!
  {
    
    #CALCULATING ERROR MEAUSRES, ordered and first selected n
    cc_er_ke <- cor(rg_cc[order(rg_cc, decreasing=TRUE)][1:te], sg_cc[order(rg_cc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    bc_er_ke <- cor(rg_bc[order(rg_bc, decreasing=TRUE)][1:te], sg_bc[order(rg_bc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    ec_er_ke <- cor(rg_ec[order(rg_ec, decreasing=TRUE)][1:te], sg_ec[order(rg_ec, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    hc_er_ke <- cor(rg_hc[order(rg_hc, decreasing=TRUE)][1:te], sg_hc[order(rg_hc, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    pr_er_ke <- cor(rg_pr[order(rg_pr, decreasing=TRUE)][1:te], sg_pr[order(rg_pr, decreasing=TRUE)][1:te]
                    , use = "pairwise.complete.obs", method = "kendall")
    
    res1 <- data.frame(
      k=k
      ,te =te
      ,ter = te / nesg2
      ,nedt = nesg2 #number of nodes of the random graph, has to be 141
      ,cc_er_ke
      ,bc_er_ke
      ,ec_er_ke
      ,hc_er_ke
      ,pr_er_ke
      ,tim=as.character(Sys.time())
    )
    
    be_sd_res<-rbind(be_sd_res,res1)
    
  }
  
  cat("k=", k, "\n")   
  
  Sys.sleep(0.01)
} 

colnames(be_sd_res)

#results important to note, that it start with node=14
be_sdt_res <-be_sd_res %>%
  group_by(ter) %>%
  summarize(med_cc_er_ke = median(cc_er_ke),
            med_bc_er_ke = median(bc_er_ke),
            med_ec_er_ke = median(ec_er_ke),
            med_hc_er_ke = median(hc_er_ke),
            med_pr_er_ke = median(pr_er_ke))

boxplot(be_sdt_res[be_sdt_res$ter==1, 2:6])
boxplot(be_sdt_res[be_sdt_res$ter<1, 2:6])
summary(be_sdt_res[,2:6], na.rm=TRUE)

baseline_error_seq_schuf <- 0 # of course +-0.01

#05.VISUALIZING the RESULTS of RAW tau-beta and baselines
#05.VISUALIZING the RESULTS of RAW tau-beta and baselines
res_fin <- res %>% 
  mutate(cc_er_ke_diff = 1-(cc_er_ke),
         bc_er_ke_diff = 1-(bc_er_ke),
         ec_er_ke_diff = 1-(ec_er_ke),
         hc_er_ke_diff = 1-(hc_er_ke),
         pr_er_ke_diff = 1-(pr_er_ke))

#first of definitive tau value for 100% of the nodes
res_ter1 <- res_fin[res_fin$ter==1,c(20:24)]
summary(res_ter1)
res1_long <-gather(res_ter1, meas, tau)

summary(res1_long$tau)

ggplot(res1_long, aes(x=reorder(meas,tau), y=tau)) +
  geom_boxplot() +
  ylim(0,0.8) +
  stat_summary(fun.y = mean, geom = "point", shape = 18, size = 3, color = "red") +
  stat_summary(fun.y = mean, geom = "text", vjust=-1, color="red",
               aes(label=sprintf("%.2f",..y..)),
               size = 2.5) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  labs(x = "Centrality measure", y = "(1-τb)", 
       title = "Centrality ranking order errors in a life-like zone [#central nodes=141 (100%)]") +
  scale_x_discrete(labels = c(cc_er_ke_diff = "Closeness", 
                              bc_er_ke_diff = "Betweenness",
                              ec_er_ke_diff = "Eigenvector",
                              hc_er_ke_diff = "Harmonic",
                              pr_er_ke_diff = "Page rank"))

#second of semi-definitive tau-b value for 10-99% of the nodes
res_fin <- res %>% 
  mutate(cc_er_ke_diff = 1-(cc_er_ke),
         bc_er_ke_diff = 1-(bc_er_ke),
         ec_er_ke_diff = 1-(ec_er_ke),
         hc_er_ke_diff = 1-(hc_er_ke),
         pr_er_ke_diff = 1-(pr_er_ke))

res_ter2 <- res_fin[res_fin$ter<1,c(7,14:18,20:24)]
res_long <- pivot_longer(data=res_ter2, cols=cc_er_ke_diff:pr_er_ke_diff)

ggplot(res_long, aes(x=ter, y=value, color=name, fill=name)) +
  stat_summary(fun.y = "median", geom = "point", shape = 22, size = 3, alpha=1) +
  labs(x = "Proportion of top nodes involved", y = "(1-τbe) [semi definitive τ-b]", 
       title = "(1-τbe) median values in life-like zone",
       fill = "Median", color = "") +
  ylim(0.0,1)+
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  scale_fill_discrete(labels = c("cc_er_ke_diff" = "Closeness c.", 
                                 "bc_er_ke_diff" = "Betweenness c.", 
                                 "ec_er_ke_diff" = "Eigenvector c.",
                                 "hc_er_ke_diff" = "Harmonic c.",
                                 "pr_er_ke_diff" = "Page rank")) +
  guides(color="none")
