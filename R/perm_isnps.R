# Synchronised permutations of a data set
#
# @param dat n x p data set matrix
#
# @return This function returns a matrix containing the partial test statistics of each pairwise comparison for each permutation applied to the data set
#
# @author Alessandra Cabassi \email{alessandra.cabassi@mail.polimi.it}
#

perm.isnps = function(dat, grp, iter, dist, load = FALSE){

    table_groups = table(grp) # groups table
    C = length(table_groups)  # number of groups

    for(i in 1:C){
        if(table_groups[i]!=table_groups[1])
            stop("All groups must be of the same size in order to apply ISNPS permutations.")
    }

    n = (dim(dat)[1])/C       # number of observations in each group
    K = C*(C-1)/2             # number of partial tests
    p = dim(dat)[2]           # number of samplings per function

    T = array(5, dim = c(ceiling(iter/100)*100+1,K)) # test statistics vector initialisation

    ### Step 2: compute test statistic

    cont_pair = 1
    for(i in 1:(C-1)){
        for(j in (i+1):C){ # for each pair of groups
            # compute test statistic for initial data
            T[1,cont_pair] = distCov(cov(dat[grp==i,],use='pairwise'),cov(dat[grp==j,],use='pairwise'),dist)
            cont_pair = cont_pair+1
        }
    }

    ### Steps 3 and 4: apply iter permutations and compute the test statistic for each permuted data set

    if(load) pb = txtProgressBar(min = 0, max = ceiling(iter/100)*100+1, style = 3) # create progress bar

    # Apply permutation to each group separately

    cont_perm = 2
    for(first_perm in 1:ceiling(iter/100)){

        perm_dat = dat
        for(i in 1:C){
            orig_grp_i = dat[grp==i,]
            perm_dat[grp==i,] = orig_grp_i[sample(n),]
        }

        # Build pseudomatrix

        X = array(0,dim = c((2*n),p,K)) # matrix 2*n X p X number of groups
        cont_pseudo = 1
        for(i in 1:(C-1)){ # for each pair of groups
            for(j in (i+1):C){
                X[,,cont_pseudo] = rbind(perm_dat[grp==i,],perm_dat[grp==j,]) # fill in the matrix
                cont_pseudo = cont_pseudo + 1
            }
        }


        # Apply 'iter' permutations to the pseudomatrix and compute test statistics

        for(bb in 1:100){
            exchanges = rbinom(n, 1, 0.5) # select permutation
            X.perm = X
            for(obs in 1:n){
                if(exchanges[obs]==0){
                    X.perm[obs,,] = X[obs,,]
                    X.perm[obs+n,,] = X[obs+n,,]
                }else{
                    X.perm[obs,,] = X[obs+n,,]
                    X.perm[obs+n,,] = X[obs,,]
                }
            }
            cont_pair = 1
            for(i in 1:(C-1)){ # for each pair of groups
                for(j in (i+1):C){ # compute test statistic for the permuted dataset
                    T[cont_perm,cont_pair] = distCov(cov(X.perm[c(1:n),,cont_pair],use = 'pairwise'),cov(X.perm[-c(1:n),,cont_pair],use = 'pairwise'),dist)
                    cont_pair = cont_pair+1
                }
            }
            cont_perm = cont_perm + 1
        } # end iter
        if(load) setTxtProgressBar(pb, cont_perm-1) # update progress bar
    }

    if(load) close(pb) # close progress bar

    return(T)
}
