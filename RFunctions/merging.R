# functions to merge the branches
# 

# auxiliary function to detect whether the leaves are associated to significant
# abundances or not
filter_abund <- function(x, thres, thres_percent)
{
    number_zeros <- length(which(x<thres))
    selected <- number_zeros <  length(x)*(thres_percent)
}

#  first step: merge brothers (if no parents of higher level is present). 
#  Iterative procedure: every time a merging is done the tree is analysed from 
#  the beginning

merge_brothers <- function(tree_old) # tree_old is a phyloseq-class experiment-level 
    # object (containing a phylogenetic tree and an otu table)
{
    nleaf <- length(phy_tree(tree_old)$tip.label)
    i <- 1
    filter_data <- apply(otu_table(tree_old), 1, function(x){filter_abund(x, 5, 0.9)})
    
    cat("merging")
    
    while(i < nleaf)
    {
        # print(i)
        
        if (!filter_data[i]) # i need to be merged
        {
            ancestor <- phy_tree(tree_old)$edge[ which(phy_tree(tree_old)$edge[,2] == i), 1]
            pair <- phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==ancestor),2]
            if(length(pair)<=1)
            {
                skip_next_level <- TRUE
                i <- i+1
            }else
            {
                if(length(which(pair>nleaf))!=0)
                {
                    skip_next_level <- TRUE
                    i <- i+1
                }else
                {
                    pair <- pair[which(pair != i)]
                    to_merge <- pair[which.min(abs(pair-i))]
                    
                    tree_new <- merge_taxa(tree_old, c(i,to_merge), archetype = 2L) #, 
                    pos_to_change <- which(phy_tree(tree_new)$tip.label == 
                                               phy_tree(tree_old)$tip.label[to_merge])
                    labels <- phy_tree(tree_new)$tip.label
                    labels[pos_to_change] <- 
                        paste(phy_tree(tree_old)$tip.label[to_merge], phy_tree(tree_old)$tip.label[i], 
                              sep = '+')
                    taxa_names(tree_new) <- c(labels)
                    otu_table(tree_new)[pos_to_change,] <- 
                        colSums(otu_table(tree_old)[c(i, to_merge),])
                       # print(length(otu_table(tree_new)[pos_to_change,]))
                        #print(length(colSums(otu_table(tree_old)[c(i, to_merge),])))
                    tree_old <- tree_new
                    # if something is merged we need to start again, but with a smaller tree
                    nleaf = nleaf - 1
                    cat('.')
                    i = 1
                    filter_data <- apply(otu_table(tree_old), 1, function(x){filter_abund(x, 5, 0.9)})
                }
            }
            
        }else
        {
            i <- i+1
        }
    }
    
    return(tree_old)
}

# then we can allow the merging with a brother also if there is a brother at an
# upper level
# 

merge_brother_with_uncle <- function(tree_old)
{
    nleaf <- length(phy_tree(tree_old)$tip.label)
    
    i <- 1
    filter_data <- apply(otu_table(tree_old), 1, function(x){filter_abund(x, 5, 0.9)})
    
    cat("merging \n")
    
    while(i <= nleaf)
    {
        if (!filter_data[i]) # i need to be merged
        {
            ancestor <- phy_tree(tree_old)$edge[ which(phy_tree(tree_old)$edge[,2] == i), 1]
            pair <- phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==ancestor),2]
            if(length(pair)<=1)
            {
                skip_next_level <- TRUE
                i <- i+1
                #print('no brothers')
            }else
            {
                if(length(which(pair>nleaf))!=0)
                {
                    pair <- pair[which(pair != i)]
                    #                 pair_small <- pair[which(pair < nleaf)]
                    #                 if(length(pair_small) == 0)
                    #                 {
                    #                     skip_next_level <- TRUE
                    #                     i <- i+1
                    #                 }else
                    #                 {
                    if(length(which(pair<nleaf))!=0)
                    {
                        pair <- pair[which(pair<nleaf)]
                    }
                    to_merge <- pair[which.min(abs(pair-i))]
                    if (to_merge<nleaf)
                    {
                        tree_new <- merge_taxa(tree_old, c(i, to_merge), archetype = 2L) #, 
                        pos_to_change <- which(phy_tree(tree_new)$tip.label == 
                                                   phy_tree(tree_old)$tip.label[to_merge])
                        labels <- phy_tree(tree_new)$tip.label
                        labels[pos_to_change] <- 
                            paste(phy_tree(tree_old)$tip.label[to_merge], phy_tree(tree_old)$tip.label[i], 
                                  sep = '+')
                        taxa_names(tree_new) <- c(labels)
                        otu_table(tree_new)[pos_to_change,] <- 
                            colSums(otu_table(tree_old)[c(i, to_merge),])

                           # print(length(otu_table(tree_new)[pos_to_change,]))
                            #print(length(colSums(otu_table(tree_old)[c(i, to_merge),])))
                        
                        tree_old <- tree_new
                        # if something is merged we need to start again, but with a smaller tree
                        nleaf = nleaf - 1
                        # print(paste("M", nleaf, sep = ' ')) 
                        cat(".")
                        i = 1
                        filter_data <- apply(otu_table(tree_old), 1, function(x){filter_abund(x, 5, 0.9)})
                    }else
                    {
                        # if to_merge is greater than nleaf => for sure the min selected
                        # is the first leaf (even if it is more farther then the others)
                        
                        # merge with the son of the brother (all the elements of pair are fathers.
                        # here I select all the sons)
                        all_the_son <- unlist(sapply(pair, function(x){
                            phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==x),2]}))
                        #  son <- phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==to_merge),2]
                        son_selected <- all_the_son[which.min(abs(all_the_son-i))] # son (!= to_merge) cannot contain i
                            
                        if( son_selected > nleaf)
                        {
                            # or with the son of the son of the brother (and stop)
                            all_the_son_son <- as.vector(sapply(all_the_son, function(x){
                                phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==x),2]}))
                            
                            son_selected <- all_the_son_son[which.min(abs(all_the_son_son-i))] # son (!= to_merge) cannot contain i
                            
                            #son <- phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==son_selected),2]
                            #son_selected <- son[which.min(abs(son-i))]
                            
                        }
                        if(son_selected>nleaf)
                        {
                            # print("Strange conf: all the brother are at most grandparents")
                            skip_next_level <- TRUE
                            i <- i+1
                        }
                        
                        tree_new <- merge_taxa(tree_old, c(i,son_selected), archetype = 2L) #, 
                        pos_to_change <- which(phy_tree(tree_new)$tip.label == 
                                                   phy_tree(tree_old)$tip.label[son_selected])
                        labels <- phy_tree(tree_new)$tip.label
                        labels[pos_to_change] <- 
                            paste(phy_tree(tree_old)$tip.label[son_selected], phy_tree(tree_old)$tip.label[i], 
                                  sep = '+')
                        taxa_names(tree_new) <- c(labels)
                        otu_table(tree_new)[pos_to_change,] <- 
                            colSums(otu_table(tree_old)[c(i, son_selected),])

                    #        print(length(otu_table(tree_new)[pos_to_change,]))
                     #       print(length(colSums(otu_table(tree_old)[c(i, son_selected),])))
                        tree_old <- tree_new
                        # if something is merged we need to start again, but with a smaller tree
                        nleaf = nleaf - 1
                        # print(paste("M", nleaf, sep = ' ')) 
                        cat('.')
                        i = 1
                        filter_data <- apply(otu_table(tree_old), 1, function(x){filter_abund(x, 5, 0.9)})
                        
                    }
                    #                 }
                }else
                {
                    pair <- pair[which(pair != i)]
                    to_merge <- which.min(abs(pair-i))
                    
                    tree_new <- merge_taxa(tree_old, c(i,to_merge), archetype = 2L) #, 
                    pos_to_change <- which(phy_tree(tree_new)$tip.label == 
                                               phy_tree(tree_old)$tip.label[to_merge])
                    labels <- phy_tree(tree_new)$tip.label
                    labels[pos_to_change] <- 
                        paste(phy_tree(tree_old)$tip.label[to_merge], phy_tree(tree_old)$tip.label[i], 
                              sep = '+')
                    taxa_names(tree_new) <- c(labels)
                    otu_table(tree_new)[pos_to_change,] <- 
                        colSums(otu_table(tree_old)[c(i, to_merge)])

                     #   print(length(otu_table(tree_new)[pos_to_change,]))
                      #  print(length(colSums(otu_table(tree_old)[c(i, to_merge),])))
                    tree_old <- tree_new
                    # if something is merged we need to start again, but with a smaller tree
                    nleaf = nleaf - 1
                   # print(paste("M", nleaf, sep = ' ')) 
                    cat('.')
                    i = 1
                    filter_data <- apply(otu_table(tree_old), 1, function(x){filter_abund(x, 5, 0.9)})
                }
            }
            
        }else
        {
            i <- i+1
        }
    }
    return(tree_old)
    
}


# merge the branches based on correlation

filter_tree_correlation <- function(tree, thres_corr, scaled = FALSE)
{
    global_tree_taxa_3 <- tree
    stop_criteria <- FALSE
    if (!scaled) {
        global_tree_taxa_5 <- phyloseq(phy_tree = phy_tree(global_tree_taxa_3), 
                                        otu_table= otu_table(t(scale(t(otu_table(global_tree_taxa_3)))), taxa_are_rows = TRUE))
        global_tree_taxa_3 <- global_tree_taxa_5
    }
    while(!stop_criteria)
    {
        edges <- phy_tree(global_tree_taxa_3)$edge
        thres_corr <- 0.7
        num_leaf <- length(phy_tree(global_tree_taxa_3)$tip.label)
        merge <- FALSE
        num_node <- num_leaf + 1
        while(!merge)
        {
            pos_merge <- edges[which(edges[,1]==num_node & edges[,2]<=num_leaf),2]
            pos_merge
            if(length(pos_merge)>=2)
            {
                
                to_analize <- combn(pos_merge, 2, simplify = FALSE)
                vector_for_merge <- rep(FALSE, num_leaf)
                sign_correlation <- rep(NA, num_leaf)
                for (i in 1:length(to_analize))
                {
                    to_analize_here <- to_analize[[i]]
                    correlation <- cor(as.vector(otu_table(global_tree_taxa_3)[to_analize_here[1],]), 
                                       as.vector(otu_table(global_tree_taxa_3)[to_analize_here[2],]))
                    if (abs(correlation) > thres_corr)
                    {
                        merge <- TRUE
                        vector_for_merge[to_analize_here] <- TRUE
                        if (correlation > 0)
                        {
                            if(length(is.na(sign_correlation[to_analize_here]))!=2)
                            {
                                if(length(which(sign_correlation[to_analize_here])=='-')!=0)
                                {
                                    stop('Positive and Negative correlation in the same branch')
                                }
                            }
                            sign_correlation[to_analize_here] <- '+'
                        }else
                        {
                            if(length(is.na(sign_correlation[to_analize_here]))!=2)
                            {
                                if(length(which(sign_correlation[to_analize_here])=='+')!=0)
                                {
                                    stop('Positive and Negative correlation in the same branch')
                                }
                            }
                            sign_correlation[to_analize_here] <- '-'
                        }
                    }
                }
                if(!merge)
                {
                    if(num_node == max(phy_tree(global_tree_taxa_3)$edge))
                    {
                        merge <- TRUE #exit because anything found
                    }else
                    {
                        num_node <- num_node + 1 
                    }
                }
                
            }else
            {
                if(num_node == max(phy_tree(global_tree_taxa_3)$edge))
                {
                    merge <- TRUE # exit because no final branches found
                }else
                {
                    num_node <- num_node + 1 
                }
            }
            
        }
        if(length(which(vector_for_merge))==0)
        {
            stop_criteria <- TRUE
        }else
        {
            global_tree_taxa_4 <- merge_taxa(global_tree_taxa_3, 
                                             which(vector_for_merge), archetype = 1L) #, 
            pos_to_change <- which(phy_tree(global_tree_taxa_4)$tip.label == 
                                       phy_tree(global_tree_taxa_3)$tip.label[which(vector_for_merge)[1]])
            if(sign_correlation[which(vector_for_merge)[1]]=='+')
            {
                labels <- phy_tree(global_tree_taxa_4)$tip.label
                labels[pos_to_change] <- 
                    paste(phy_tree(global_tree_taxa_3)$tip.label[which(vector_for_merge)], 
                          collapse = '&')
                taxa_names(global_tree_taxa_4) <- c(labels)
                
                vec_here <- colMeans(otu_table(global_tree_taxa_3)[which(vector_for_merge),])
                vec_here_scaled <- (vec_here - mean(vec_here))/sd(vec_here)
                otu_table(global_tree_taxa_4)[pos_to_change,] <- vec_here_scaled
                   
            
              #      print(length(colMeans(otu_table(global_tree_taxa_3)[which(vector_for_merge),])))
               #     print(length(otu_table(global_tree_taxa_4)[pos_to_change,]))
                    
            }
            if(sign_correlation[which(vector_for_merge)[1]]=='-')
            {
                print('Negative correlation')
                if(length(which(vector_for_merge))!=2)
                {
                    stop('more than 2 elements with negative correlation')
                }
                labels <- phy_tree(global_tree_taxa_4)$tip.label
                labels[pos_to_change] <- 
                    paste(phy_tree(global_tree_taxa_3)$tip.label[which(vector_for_merge)], 
                          collapse = '%')
                taxa_names(global_tree_taxa_4) <- c(labels)
                
                vec_here <- (otu_table(global_tree_taxa_3)[which(vector_for_merge)[1]]
                     - otu_table(global_tree_taxa_3)[which(vector_for_merge)[2]])/2
                vec_here_scaled <- (vec_here - mean(vec_here))/sd(vec_here)
                otu_table(global_tree_taxa_4)[pos_to_change,] <- vec_here_scaled
                    
                
            }
          #  global_tree_taxa_5 <- transform_sample_counts(global_tree_taxa_4, 
          #                                                function(x) {(x-mean(x)/sd(x))})
          
         # global_tree_taxa_5 <- phyloseq(phy_tree = phy_tree(global_tree_taxa_3), 
          #                                    otu_table= otu_table(t(scale(t(otu_table(global_tree_taxa_3)))), taxa_are_rows = TRUE))
          global_tree_taxa_3 <- global_tree_taxa_4
        }
    }
    return(global_tree_taxa_3)
}



# remove from a tree the sons of a set of internla nodes
#

remove_son <- function(tree_old, tab_out)
{
    names_removed <- NULL
    cat("removing \n")
    
    
    bad_start <- which(!(rownames(tab_out) %in% tree_old$tip.label))
    no_leaves <- rownames(tab_out[bad_start, ])
    middle_knots <- tree_old$node.label %in% no_leaves
    print(table(middle_knots))
    middle_knots_to_remove <- which(middle_knots)
    names_knots_to_remove <- tree_old$node.label[middle_knots]

    while (length(middle_knots_to_remove) >= 1)
    {
      #  print(table(middle_knots))
        cat(".")
        
        pos_knot_to_remove <- which(tree_old$node.label==names_knots_to_remove[1]) + 
            length(tree_old$tip.label)
       # print(pos_knot_to_remove)
        son <- tree_old$edge[which(tree_old$edge[,1]==pos_knot_to_remove),2]
     #   print(son)
      #  print(length(tree_old$tip.label))
        if(length(son)>1)
        {
            names_son <- tree_old$tip.label[son]
            
            tree_new <- merge_taxa(tree_old, son, archetype = 1L) #,
            #tree_new <- drop.tip(tree_old, son, trim.internal = FALSE)
            pos_to_change <- which(phy_tree(tree_new)$tip.label == names_son[1])
            labels <- phy_tree(tree_new)$tip.label
            labels[pos_to_change] <- names_knots_to_remove[1]
            taxa_names(tree_new) <- labels
        }
        if(length(son)==1)
        {
            tree_new <- tree_old
            pos_to_change <- which(phy_tree(tree_new)$tip.label == tree_old$tip.label[son])
            labels <- phy_tree(tree_new)$tip.label
            labels[pos_to_change] <- names_knots_to_remove[1]
            taxa_names(tree_new) <- labels
        }
        if(length(son)==0)
        {
            print('no son')
        }
      
        tree_old <- collapse.singles(tree_new)
        names_removed <- c(names_removed, names_son)
        
        # print(length(middle_knots_to_remove))
        bad_start <- which(!(rownames(tab_out) %in% tree_old$tip.label))
        no_leaves <- rownames(tab_out[bad_start, ])
        middle_knots <- tree_old$node.label %in% no_leaves
        table(middle_knots)
        middle_knots_to_remove <- which(middle_knots)
        names_knots_to_remove <- tree_old$node.label[middle_knots]
        length(names_knots_to_remove)
    }

    return(list(tree = tree_old, names_removed = names_removed))
}

#tree_old <- global_tree_taxa_3

detect_sons <- function(x, tree_old = tree_old)
{
    phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==x),2]
}

# tree_old <- global_tree_taxa_3
cut_phyla <- function(tree_old)
{
    detect_sons <- function(x, tree_old = tree_old)
    {
        phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,1]==x),2]
    }
    
    ended_phyla <- FALSE
    while ( !ended_phyla )
    {
        # define the root
        ances <- phy_tree(tree_old)$edge[which(phy_tree(tree_old)$edge[,2]==1),1]
        while (length(ances)!=0)
        {
           # print(ances)
            ances_old <- ances
            ances <- phy_tree(tree_old)$edge[
                which(phy_tree(tree_old)$edge[,2]==ances_old),1]
        }
       # print(ances)
        
         # print(ances_old)
        
        root <- ances_old
        number_lab <- length(phy_tree(tree_old)$tip.label)
       # print(number_lab)
        name_root <- phy_tree(tree_old)$node.label[root - number_lab]
       # print(name_root)
        
        root_son_2 <- phy_tree(tree_old)$edge[
            which(phy_tree(tree_old)$edge[,1]==root),2]
        
        names_phyla <- unlist(sapply(root_son_2, function(x){
            if (x > number_lab)
            {
                return( phy_tree(tree_old)$node.label[x - number_lab])
            }else
            {
                return(  phy_tree(tree_old)$tip.label[x] )
            }
        }))
        
        t <- 1
        r <- root_son_2[1]
        number_final_phyla <- NULL
        
        if (r < number_lab)
        {
            number_final_phyla <- 1
        }else
        {
            number_final_phyla <- 0
        }
        
        while( (r < number_lab) & (t < length(names_phyla)))
        {
            r <-  root_son_2[t+1]
            t <- t+1
            number_final_phyla <- number_final_phyla + 1
        }
        
        if (number_final_phyla == length( names_phyla))
        {
            ended_phyla <-  TRUE
        }
        
        if (!ended_phyla)
        {
            to_be_merged <- phy_tree(tree_old)$edge[
                which(phy_tree(tree_old)$edge[,1]==r),2]
            while(length(which((to_be_merged>number_lab))) != 0) # while some son is not a leaf
            {
              #  cat(".")
                to_be_merged_to_be_added <- sapply(to_be_merged, detect_sons, tree_old = tree_old)
             #   print(to_be_merged_to_be_added)
                to_be_merged_new <- sapply( 1:length(to_be_merged_to_be_added), 
                                            function(x){
                                                if (length(to_be_merged_to_be_added[[x]])!=0)
                                                {
                                                    return(to_be_merged_to_be_added[[x]])
                                                }else
                                                {
                                                    return(to_be_merged[[x]])
                                                }})
                to_be_merged <- unlist(to_be_merged_new)
               # print(to_be_merged)
            }
           # print(to_be_merged)
            if (length(to_be_merged)!=0) 
            {
                tree_new <- merge_taxa(tree_old, to_be_merged, archetype = 1L )
                
                pos_to_change <- which(phy_tree(tree_new)$tip.label == 
                                           phy_tree(tree_old)$tip.label[to_be_merged[1]])
                
                otu_table(tree_new)[pos_to_change,] <- 
                    colSums(otu_table(tree_old)[to_be_merged,])
                
                labels <- phy_tree(tree_new)$tip.label
                labels[pos_to_change] <- names_phyla[t]
                taxa_names(tree_new) <- labels
            }
            tree_old <- tree_new
        }
       
    }
    
    return(tree_old)
}

