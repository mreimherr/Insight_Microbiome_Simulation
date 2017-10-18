# This code goes through all of the steps to 
# produce the final table of abundances used in FLAME
# Final values are at very end, where a latex table of labels
# and a data matrix are produced

# Various plots that can be produced are included, but commented out

rm(list = ls())
library(ctv) 
library(ape)
library(phyloseq)
library(knitr)


new_abundaces <- read.csv('data/abundances_Mar2017.csv')
names_ab <- colnames(new_abundaces)[-(1:3)]
colnames(new_abundaces)[1:3]


new_abundaces_Cstool <- new_abundaces[which((new_abundaces$type == 'stool' | new_abundaces$type == 'Stool') & new_abundaces$who=='C1'),]


# 115 e 269 don't have growth measurement => removed
new_abundaces_Cstool$ID[c(12,120)]

new_abundaces_Cstool <- new_abundaces_Cstool[-c(12,120),]
dim(new_abundaces_Cstool)
list_IDs_Cstool <- new_abundaces_Cstool$ID

new_abundaces_sort_Cstool <- new_abundaces_Cstool[with(new_abundaces_Cstool, 
                                                         order(ID)), ]

new_abundaces_Cstool <- new_abundaces_sort_Cstool
dim(new_abundaces_Cstool)
list_IDs_Cstool <- new_abundaces_Cstool$ID


library(phytools)
tree <- phytools::read.newick(file="data/NewickTree_15Sept2016.txt")

tree<-collapse.singles(tree)

# # pdf('new_tree_15Sept.pdf', height = 200, width = 20)
# # par(mar=c(0,0,0,0))
# # plot(tree)
# # dev.off()

head(tree$edge)
head(tree$tip.label)



# Child Stool ------------------------------------------------------------

# Analysis of the Child Stool Abundances

tab_out <- otu_table(t(new_abundaces_Cstool[,-(1:3)]), TRUE)

length(which(tree$tip.label %in% rownames(tab_out)))
# number of bacteria of the dataset present also in the tree
length(which(rownames(tab_out) %in%  tree$tip.label ))

library(stringr)

# modify names of the tree to adapt to the names of abundances
# tip
new_treelbel <- str_replace_all(tree$tip.label, "[']", "")
tree$tip.label <- new_treelbel
# node
new_nodelabel <- str_replace_all(tree$node.label, "[']", "")
tree$node.label <- new_nodelabel

start_str <- rownames(tab_out) # names of abundances (also str_2, str_3... 
# are names of abundances after changes to adapt to the tree)
length(which(tree$tip.label %in% start_str))

# modify names of the abundaces to adapt to the names of the tree
str_2 <- str_replace_all(start_str, "[.]", "-")
length(which(tree$tip.label %in% str_2))

str_3 <- str_2
bad_start <- which(!(str_3 %in% tree$tip.label))
str_3[bad_start] <- str_replace_all(str_3[bad_start], "[-]", "_")
# # now the bacteria present in the two samples are
length(which(tree$tip.label %in% str_3))

rownames(tab_out) <- str_3

source("RFunctions/merging.R")
str_4 <- str_3
bad_start <- which(!(str_4 %in% tree$tip.label))
no_leaves <- rownames(tab_out[bad_start, ])


middle_knots <- tree$node.label %in% no_leaves
table(middle_knots) # 88 abundances are not leaves but middle knots of the tree
head(tree$node.label[middle_knots])

middle_knots_to_remove <- which(middle_knots)

tree_new_here <- remove_son(tree, tab_out)
removed_names <- tree_new_here$names_removed
tree_new <- tree_new_here$tree
# # 
# # pdf('new_tree_15Sept_selected.pdf', height = 200, width = 20)
# # par(mar=c(0,0,0,0))
# # plot(tree_new)
# # dev.off()
# # 
# # 
# # pdf('new_tree_15Sept_names.pdf', height = 200, width = 20)
# # par(mar=c(0,0,0,0))
# # plot(tree, show.node.label=TRUE)
# # dev.off()

length(which(removed_names %in% start_str))
# initial defintion of the tree: leaves that are abundances
length(which(tree$tip.label %in% start_str))
# initial defintion of the tree: leaves that are abundances after string modifications
length(which(tree$tip.label %in% str_3))
# initial defintion of the tree: middle knots that are abundances
length(which(tree$node.label %in% str_3))
# after the merging: leaves that are abundances
length(which(tree_new$tip.label %in% str_3))
length(which(tree_new$node.label %in% str_3)) # ok: all the internal nodes are leaves

# analysis of names with double names
names_tips <- tree_new$tip.label
head(names_tips)
first_parts <- sapply(names_tips, function(x){strsplit(x, "_", 1)})
# position of elements with two names
# and substitution of the names of abundances with
# _ in the name of the tree, but present only one time
twice <- sapply(first_parts, function(x){length(x)>1})
first_elements <- sapply(first_parts, function(x){x[1]})
number_elements <- sapply(first_elements, function(x){ 
    vec_presence <- sapply(tree_new$tip.label, 
                           function(y){length(grep(x,y))})
    return(sum(vec_presence))})
names_tips_2 <- ifelse(number_elements==1, first_elements, names_tips )

length(which(names_tips %in% start_str))
length(which(names_tips %in% str_3))
length(which(names_tips_2 %in% str_3))

length(start_str[which(!(str_3 %in% names_tips_2))])
start_str[which(!(str_3 %in% names_tips_2))]

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Blochmannia",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Blochmannia"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Methylacidiphilum",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Methylacidiphilum"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Portiera",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Portiera"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Neoehrlichia",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Neoehrlichia"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Koribacter",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Koribacter"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Aquiluna",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Aquiluna"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Glomeribacter",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Glomeribacter"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Microthrix",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Microthrix"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Regiella",x))})==1)
names_tips_2[pos_blo] <- "Candidatus_Regiella"

length(start_str[which(!(str_3 %in% names_tips_2))])
start_str[which(!(str_3 %in% names_tips_2))]

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("1-68",x))})==1)
names_tips_2[pos_blo] <- "X1_68"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("02d06",x))})==1)
names_tips_2[pos_blo] <- "X02d06"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("31d11",x))})==1)
names_tips_2[pos_blo] <- "X31d11"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("7N15",x))})==1)
names_tips_2[pos_blo] <- "X5_7N15"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("A55",x))})==1)
names_tips_2[pos_blo] <- "A55_D21"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("aeronauticum",x))})==1)
names_tips_2[pos_blo] <- "Desulfotomaculum"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Fibrobacter",x))})==1)
pos <- which(sapply(start_str, function(x){length(grep("Fibrobacter",x))})==1)

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("succinogenes",x))})==1)
names_tips_2[pos_blo] <- "Fibrobacter" # one of the two Fibrobacter has already a name _2

length(start_str[which(!(str_3 %in% names_tips_2))])
start_str[which(!(str_3 %in% names_tips_2))]

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Geosporobacter",x))})==1)
names_tips_2[pos_blo] <- "Geosporobacter_Thermotalea"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("L7A",x))})==1)
names_tips_2[pos_blo] <- "L7A_E11"


pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Natronincola",x))})==1)
names_tips_2[pos_blo] <- "Natronincola_Anaerovirgula"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Pseudoramibacter",x))})==1)
names_tips_2[pos_blo] <- "Pseudoramibacter_Eubacterium"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Tindallia",x))})==1)
names_tips_2[pos_blo] <- "Tindallia_Anoxynatronum"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("Tissierella",x))})==1)
names_tips_2[pos_blo] <- "Tissierella_Soehngenia"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("WAL",x))})==1)
names_tips_2[pos_blo] <- "WAL_1855D"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("hetero",x))})==1)
names_tips_2[pos_blo] <- "heteroC45_4W"

pos_blo <- which(sapply(names_tips_2, function(x){length(grep("aceti",x))})==1)
names_tips_2[pos_blo] <- "Acetobacter"


pos <- which(sapply(start_str, function(x){length(grep("Prevotella",x))})==1)
pos
pos_blo <- which(sapply(names_tips_2, function(x){length(grep("tannerae",x))})==1)
names_tips_2[pos_blo] <- "[Prevotella]"

length(start_str[which(!(str_3 %in% names_tips_2))])

tree_new$tip.label <- names_tips_2
length(which(duplicated(names_tips_2)))

# two double matches. 
# 1) MSBL3 is present twice in the tree
# duplicated(new_tree$tip.label) => removed
pos_taxa <- sapply(tree_new$tip.label, function(x){length(grep("MSBL3",x))})!=1
table(pos_taxa)
which(pos_taxa == TRUE)
new_tree_2 <- prune_taxa(pos_taxa, tree_new)

# MSBL3 is removed also from the table
pos_blo <- which(sapply(rownames(tab_out), function(x){length(grep("MSBL3",x))})==1)
tab_out2 <- tab_out[- pos_blo,]

# 2) 4-29 is present three times in the tree
# removed from the tree
pos_taxa <- sapply(new_tree_2$tip.label, function(x){length(grep("4-29",x))})!=1
table(pos_taxa)
new_tree_3 <- prune_taxa(pos_taxa, new_tree_2)

# and from the table
pos_blo <- which(sapply(rownames(tab_out2), function(x){length(grep("4_29",x))})==1)
tab_out3 <- tab_out2[- pos_blo,]

# 3) Prevotella is present twice in the tree
pos_taxa <- sapply(new_tree_3$tip.label, function(x){length(grep("Prevotella",x))})!=1
table(pos_taxa)
new_tree_4 <- prune_taxa(pos_taxa, new_tree_3)

# and from the table
pos_blo <- which(sapply(rownames(tab_out3), function(x){length(grep("Prevotella",x))})==1)
tab_out4 <- tab_out3[- pos_blo,]



# 4) Eubacterium is present twice in the tree
pos_taxa <- sapply(new_tree_4$tip.label, function(x){length(grep("Eubacterium",x))})!=1
table(pos_taxa)
new_tree_5 <- prune_taxa(pos_taxa, new_tree_4)

# and from the table
pos_blo <- which(sapply(rownames(tab_out4), function(x){length(grep("Eubacterium",x))})==1)
tab_out5 <- tab_out4[- pos_blo,]


# 5) Ruminococcus is present twice in the tree
pos_taxa <- sapply(new_tree_4$tip.label, function(x){length(grep("Ruminococcus",x))})!=1
table(pos_taxa)
new_tree_5 <- prune_taxa(pos_taxa, new_tree_4)

# and from the table
pos_blo <- which(sapply(rownames(tab_out4), function(x){length(grep("Ruminococcus",x))})==1)
tab_out5 <- tab_out4[- pos_blo,]



dim(tab_out5)
new_tree_5

# # # 
# #  pdf('new_tree_15Sept_final_selection.pdf', height = 200, width = 20)
# #  par(mar=c(0,0,0,0))
# #  plot(new_tree_3)
# #  dev.off()



global_tree_taxa <- phyloseq(phy_tree = new_tree_5, otu_table = tab_out5)

global_tree_taxa

# # 
# #  pdf('new_tree_15Sept_subset_taxa.pdf', height = 200, width = 20)
# #  par(mar=c(0,0,0,0))
# #  plot(phy_tree(global_tree_taxa))
# #  dev.off()


filter_abund <- function(x, thres, thres_percent)
{
    number_zeros <- length(which(x<thres))
    selected <- number_zeros <  length(x)*(thres_percent)
}


global_step_1 <- merge_brothers(global_tree_taxa)
global_step_1

filter_data <- apply(otu_table(global_step_1), 1, function(x){filter_abund(x, 5, 0.9)})
table(filter_data)

# second step

global_step_2 <- merge_brother_with_uncle(global_step_1)
global_step_2

filter_data <- apply(otu_table(global_step_2), 1, function(x){filter_abund(x, 5, 0.9)})
table(filter_data)
# 
# par(mar=c(0,0,0,0))
# plot(phy_tree(global_step_2), x.lim = 180, cex = 0.5, tip.color = ifelse(filter_data, 'darkgreen', 'red3') )

global_tree_taxa_3 <- phyloseq(phy_tree = phy_tree(global_step_2), 
                               otu_table = otu_table(t(scale(t(otu_table(global_step_2)))), taxa_are_rows = TRUE))


tab_abund <- otu_table(global_tree_taxa_3)
otu_tab_gt3 <- tab_abund
rownames(otu_tab_gt3) <- 1:dim(tab_abund)[1]
phy_tree_gt3 <- phy_tree(global_tree_taxa_3)
phy_tree_gt3$tip.label <- 1:dim(tab_abund)[1]
global_tree_taxa_3_to_plot <- phyloseq(otu_table = otu_tab_gt3, phy_tree= phy_tree_gt3)

cor_matrix <- cor(t(otu_table(global_tree_taxa_3_to_plot)))

sig_mat <- ifelse(abs(cor_matrix)>0.7, 1, 0)
diag(sig_mat) <- 0



# # 
# # pdf('original_correlation_CB_2.pdf', height = 15, width = 14)
# # corrplot(cor_matrix, type="upper", order = "original", 
# #          tl.cex = 0.7, tl.col = 1, p.mat = sig_mat,  mar = c(0, 0, 1, 0), cl.cex = 2,
# #          sig.level = 0.1, pch.cex = 0.3, pch.col = 'yellow')
# # dev.off()


library(corrplot)
# corrplot(cor_matrix, type="upper", order = "original", 
# tl.cex = 0.6, tl.col = 1, p.mat = sig_mat, tl.pos='n',
# sig.level = 0.1, pch.cex = 0.2, pch.col = 'yellow')

start <- proc.time()
global_tree_taxa_4 <- filter_tree_correlation(global_tree_taxa_3,  0.7)
duration <- proc.time()-start
global_tree_taxa_4
duration

# par(mar=c(0,0,0,0))
# plot(phy_tree(global_tree_taxa_4))

# tab_abund <- otu_table(global_tree_taxa_4)
# cor_matrix <- cor(t(tab_abund))

# colnames(cor_matrix) <- 1:dim(cor_matrix)[2] 
# rownames(cor_matrix) <- 1:dim(cor_matrix)[2]


# sig_mat <- ifelse(abs(cor_matrix)>0.7, 1, 0)
# diag(sig_mat) <- 0

# corrplot(cor_matrix, type="upper", order = "original", 
# tl.cex = 0.6, tl.col = 1, p.mat = sig_mat, 
# sig.level = 0.1, pch.cex = 0.3, pch.col = 'yellow')


# start <- proc.time()
# global_tree_taxa_4_to_plot <- filter_tree_correlation(global_tree_taxa_3_to_plot,  0.7)
# duration <- proc.time()-start

# tab_abund <- otu_table(global_tree_taxa_4_to_plot)
# cor_matrix <- cor(t(tab_abund))

# sig_mat <- ifelse(abs(cor_matrix)>0.7, 1, 0)
# diag(sig_mat) <- 0
# # 
# # pdf('final_correlation_CB_2.pdf', height = 15, width = 13)
# # corrplot(cor_matrix, type="upper", order = "original", 
# #          tl.cex = 0.7, tl.col = 1, p.mat = sig_mat,  mar = c(0, 0, 1, 0), cl.cex = 2,
# #          sig.level = 0.1, pch.cex = 0.3, pch.col = 'yellow')
# # dev.off()
# # 
# 
# tab_abund <- otu_table(global_tree_taxa_4_to_plot)
# otu_tab_gt5 <- tab_abund
# 
# has_and <- grepl("&", phy_tree(global_tree_taxa_4_to_plot)$tip.label)
# name_and <- rep("", length(has_and))
# name_and[has_and] <- paste("G", 1:sum(has_and), " : ", sep = '')
# 
# names <- paste(name_and, phy_tree(global_tree_taxa_4_to_plot)$tip.label)
# 
# rownames(otu_tab_gt5) <- names
# phy_tree_gt5 <- phy_tree(global_tree_taxa_4_to_plot)
# phy_tree_gt5$tip.label <- names
# global_tree_taxa_5_to_plot <- phyloseq(otu_table = otu_tab_gt5, phy_tree= phy_tree_gt5)

# # 
# # pdf('Final_tree_CB_2.pdf', height = 15, width = 4)
# # par(mar= c(0,0,0,0))
# # plot(phy_tree(global_tree_taxa_5_to_plot), cex = 0.8)
# # dev.off()
# # 
# # 
# # pdf('Final_tree_round_CB_2.pdf', height = 10, width = 10)
# # par(mar= c(0,0,0,0), oma= c(0,0,0,0))
# # plot(phy_tree(global_tree_taxa_5_to_plot),type = "fan", cex = 0.8)
# # dev.off()

global_tree_taxa_4
al <- cbind(1:75, as.vector(phy_tree(global_tree_taxa_4)$tip.label))
head(al)



# Names for abundance groups are saved in a latex table
library(xtable)
x <- xtable(al)
print(xtable(al), type="latex", file="data/output_CS_2_March2017.tex", include.rownames=FALSE, include.colnames=FALSE)


# Table of combined abundences
# This table is scaled so that the columns have mean zero and var 1 
# before applying FLAME
tab_abund <- t(otu_table(global_tree_taxa_4))
X <- scale(tab_abund)

