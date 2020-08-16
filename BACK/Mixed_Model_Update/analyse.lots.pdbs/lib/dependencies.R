#library(bio3d)
#library(veriNA3d)

#pdblist <- dir("pdbs/", full.names=T)
#pdbfile <- pdblist[1]
#system.time(measureNuc(read.pdb(pdbfile)))
#source("measureNuc_light.R")

## Wrapper to choose between lapply and mclapply accordingly
#library(parallel)
xlapply <-
function(X, FUN, ..., mc.cores=1, mc.preschedule=FALSE) {

    if (mc.cores > detectCores()) {
        warning("The machine does not have ", mc.cores, " cores. Using ",
                    max(1, detectCores() - 1), ".")
        mc.cores <- max(1, detectCores() - 1)
    }

    if (mc.cores > 1) {
        mclapply(X=X, FUN=FUN, ...=...,
                    mc.cores=mc.cores, mc.preschedule=mc.preschedule)
    } else {
        lapply(X=X, FUN=FUN, ...=...)
    }
}

torsionals <- as.data.frame(matrix(c(
                    "pre_O3'",  "P",     "O5'",      "C5'",         "alpha",
                    "P",        "O5'",   "C5'",      "C4'",         "beta",
                    "O5'",      "C5'",   "C4'",      "C3'",         "gamma",
                    "C5'",      "C4'",   "C3'",      "O3'",         "delta",
                    "C4'",      "C3'",   "O3'",      "post_P",      "epsilon",
                    "C3'",      "O3'",   "post_P",   "post_O5'",    "zeta",
                    "O4'",      "C1'",   "N_base",   "C_base",      "chi"
                    ), ncol=5, byrow=TRUE),
                stringsAsFactors=FALSE)
colnames(torsionals) <- c("atomA", "atomB", "atomC", "atomD", "labels")

#system.time(measureNuc(read.pdb(pdbfile), angles=NULL, distances=NULL, Dp=NULL, pucker=F, torsionals=torsionals))

measureNuc_light <- function(pdb, torsionals=torsionals, refatm="C4'") {
    pdb <- veriNA3d:::.perfect_input_format(pdb)
    chain <- as.character(unique(pdb$atom$chain))
    combinations <- expand.grid(model=1, chain=chain, stringsAsFactors=FALSE)

    ntinfo <- lapply(1:nrow(combinations), function(i, torsionals) {

        chain <- combinations$chain[i]
        ## Selection of Chain of interest ------------------------------------
        selection <- atom.select(pdb, chain=chain)

        ## pdb contains the PDB object ONLY with the selected model and chain
        pdb1 <- trim(pdb, selection)

        ridlist <- pdb1$atom$resid[which(pdb1$atom$elety == c(refatm))]
        reslist <- pdb1$atom$resno[which(pdb1$atom$elety == c(refatm))]
        inslist <- pdb1$atom$insert[which(pdb1$atom$elety == c(refatm))]

        total <- length(reslist)

        indices <- seq_len(total)

        ## Call to do the maths for the given chain ------------------------------
        ntinfo <- lapply(indices,
                            FUN=new_measure,
                            reslist=reslist,
                            inslist=inslist,
                            ridlist=ridlist,
                            pdb=pdb1,
                            torsions=torsionals)

        ntinfo <- do.call(rbind, ntinfo)
        return(ntinfo)
    }, torsionals=torsionals)
    ntinfo <- do.call(rbind, ntinfo)
    return(ntinfo)
}

new_measure <-
function(index, reslist, inslist, ridlist, pdb, torsions) {
    ## Extract info about residue number, insert records and 
    ## identifier (A, G, C, U, DA, DG, DC, DT...) ----------------------------
    resno <- reslist[index]
    insert <- inslist[index]
    resid <- ridlist[index]

    ## Extract info about residue number and insert records for the previous 
    ## and following nucleotides for inter-nucleotide measures (e.g.eta-theta)
    if (index == 1) {
        preresno <- ""
        postresno <- reslist[index + 1]
        preinsert <- ""
        postinsert <- inslist[index + 1]
    }else if (index == length(reslist)) {
        preresno <- reslist[index-1]
        postresno <- ""
        preinsert <- inslist[index-1]
        postinsert <- ""
    } else {
        preresno <- reslist[index-1]
        postresno <- reslist[index + 1]
        preinsert <- inslist[index-1]
        postinsert <- inslist[index + 1]
    }

    ## Find base type --------------------------------------------------------
    if (resid %in% c("A", "G", "DA", "DG")) {
        base_type <- "pu"
    } else if (resid %in% c("C", "U", "DC", "DT")) {
        base_type <- "py"
    } else {
        ## For modified bases check the existence of the atom N9
        if (any(pdb$atom$resno == resno &
                pdb$atom$insert == insert &
                pdb$atom$elety == "N9")) {
            base_type <- "pu"
        } else {
            base_type <- "py"
        }
    }

    ## Since the input may contain some atom names called N_base and C_base,
    ## they will be replaced by N9/N1 and C4/C2 in function of base type -----
    if (base_type == "pu") {
        N_base <- "N9"
        C_base <- "C4"
    } else if (base_type == "py") {
        N_base <- "N1"
        C_base <- "C2"
    }

    torsions[torsions == "N_base"] <- N_base
    torsions[torsions == "C_base"] <- C_base

    toratoms <- unique(unlist(
                            torsions[, grep("atom", names(torsions))],
                            use.names=FALSE))

    ## Generate vector with all necessary atom names that will be used for the
    ## calculations ----------------------------------------------------------
    atomlist <- sort(unique(toratoms))

    ## Generate vector with the name of the objects that will contain the atom
    ## selection -------------------------------------------------------------
    atomlist_sel <- paste( gsub("\'", "p", atomlist), "_sel", sep="")

    ## Generate vector with real elety names (remove prefix "pre" and "post")
    atomelety <- lapply(strsplit(atomlist, split="_"),
                        function(x) {
                            return(x[length(x)])
                        })

    ## Generate vector with apropiate strings to call the correct object, 
    ## in case there are atoms of the previous or following nucleotides to be
    ## used ------------------------------------------------------------------
    prefix <- vector("character", length(atomlist))
    prefix[grep("pre", atomlist)] <- "pre"
    prefix[grep("post", atomlist)] <- "post"

    ## Generate vectors indicating which resno and insert objects should be 
    ## used ("", "pre" or "post" nucleotide records) -------------------------
    resnocall <- paste( prefix, "resno", sep="")
    insertcall <- paste( prefix, "insert", sep="")

    ## Use the previous vectors to generate as many objects as atoms, 
    ## containing the selection, the atom and xyz indices, necessary to
    ## compute distances, angles and torsionals ------------------------------
    invisible(mapply(FUN=select_many,
                        out_object=atomlist_sel,
                        elety=atomelety,
                        resno=resnocall,
                        insert=insertcall,
                        MoreArgs=list(pdb=pdb)))

    torsions$atomA <- paste(gsub("\'", "p", torsions$atomA), "_sel", sep="")
    torsions$atomB <- paste(gsub("\'", "p", torsions$atomB), "_sel", sep="")
    torsions$atomC <- paste(gsub("\'", "p", torsions$atomC), "_sel", sep="")
    torsions$atomD <- paste(gsub("\'", "p", torsions$atomD), "_sel", sep="")

    torsions_list <- paste(torsions$labels, sep ="")
    torsions_sel <- c(paste(torsions_list, "_sel", sep=""))

    invisible(apply(cbind(  torsions_sel,
                            torsions$atomA,
                            torsions$atomB,
                            torsions$atomC,
                            torsions$atomD),
                    MARGIN=1, FUN=append_selections))

    invisible(mapply(FUN=torsions_many,
                            out_object=torsions_list,
                            selection=torsions_sel,
                            MoreArgs=list(pdb=pdb)))

    toshift <- torsions_list
    invisible(lapply(toshift,
                        FUN=function(tor) {
                            assign(tor, value=shift360(get(tor)),
                                    envir=parent.frame(n=2))
                        }))
    names(torsions_list) <- torsions$labels
    torsions_list <- lapply(torsions_list, function(tor) get(tor))

    return(data.frame(torsions_list, row.names=NULL))
}

select_many <- function(out_object, elety, pdb, resno, insert) {
    resno <- get(resno, envir=parent.frame(n=2))
    insert <- get(insert, envir=parent.frame(n=2))
    assign(out_object,
            value=atom.select(pdb,
                                eleno=pdb$atom[which(
                                                pdb$atom$resno == resno &
                                                pdb$atom$insert == insert &
                                                pdb$atom$elety == elety),
                                                "eleno"],
                                verbose=FALSE),
            envir=parent.frame(n=2))
}

torsions_many <- function(out_object, selection, pdb) {
    sel <- get(selection, envir=parent.frame(n=2))
    if (length(sel) == 12) {
        assign(out_object, round(torsion.xyz(pdb$xyz[sel]), 3),
        envir=parent.frame(n=2))
    } else {
        assign(out_object, NULL, envir=parent.frame(n=2))
    }
}

append_selections <- function(selections) {
    output <- invisible(lapply(seq(2, length(selections), 1), FUN=function(i) {
        return(get(selections[i], envir=parent.frame(n=4))$xyz)
    }))
    assign(selections[1], value=unlist(output), envir=parent.frame(n=2))
}
## Function to shift 360 degrees torsion angles
shift360 <-
function(tor) {
    if (!is.null(tor) && !is.na(tor) && length(tor) > 0) {
        if (tor < 0) {
            tor_shifted <- tor + 360
        } else {
            tor_shifted <- tor
        }
        return(tor_shifted)

    } else {
        return(NA)
    }
}

scoring_bayes <- function(ntinfo_disc= NULL, bd = bayesian_data, ntinfo = NULL){
    
    #windows <- bd$windows
    alpha <- bd$alpha
    beta <- bd$beta
    gamma <- bd$gamma
    delta <- bd$delta
    epsilon <- bd$epsilon
    zeta <- bd$zeta
    chi <- bd$chi
    alpha_plus <-bd$alpha_plus
    beta_plus <- bd$beta_plus
    gamma_plus <- bd$gamma_plus
    delta_plus <- bd$delta_plus
    chi_plus <- bd$chi_plus
    beta_minus <- bd$beta_minus
    gamma_minus <- bd$gamma_minus
    delta_minus <- bd$delta_minus
    epsilon_minus <- bd$epsilon_minus
    zeta_minus <-bd$zeta_minus
    chi_minus <- bd$chi_minus
    
    #eval <- minipack::discretize(ntinfo_m = evaluate, windows = windows)
    eval <- ntinfo_disc
    eval2 <- ntinfo
    
    bd <- bayesian_data
    
    ind <-which(eval %in% 0)
    if(length(ind) > 0){
        eval[ind] <- 1
    }
    
    percent <- c()
    
    alpha_amb <- c(nulling(dnorms_calc("beta",eval[,"beta"],"alpha", eval2[,"alpha"])),
                   nulling(dnorms_calc("gamma",eval[,"gamma"],"alpha",eval2[,"alpha"])),
                   nulling(dnorms_calc("delta_minus",eval[,"delta_minus"],"alpha",eval2[,"alpha"])),
                   nulling(dnorms_calc("delta",eval[,"delta"],"alpha",eval2[,"alpha"])),
                   nulling(dnorms_calc("epsilon_minus",eval[,"epsilon_minus"],"alpha",eval2[,"alpha"])),
                   nulling(dnorms_calc("zeta_minus",eval[,"zeta_minus"],"alpha",eval2[,"alpha"])),
                   nulling(dnorms_calc("chi",eval[,"chi"],"alpha",eval2[,"alpha"])))
    
    beta_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"beta",eval2[,"beta"])),
                  nulling(dnorms_calc("gamma",eval[,"gamma"],"beta",eval2[,"beta"])),
                  nulling(dnorms_calc("delta_minus",eval[,"delta_minus"],"beta",eval2[,"beta"])),
                  nulling(dnorms_calc("delta",eval[,"delta"],"beta",eval2[,"beta"])),
                  nulling(dnorms_calc("epsilon_minus",eval[,"epsilon_minus"],"beta",eval2[,"beta"])),
                  nulling(dnorms_calc("zeta_minus",eval[,"zeta_minus"],"beta",eval2[,"beta"])),
                  nulling(dnorms_calc("chi",eval[,"chi"],"beta",eval2[,"beta"])))
    
    gamma_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"gamma",eval2[,"gamma"])),
                   nulling(dnorms_calc("beta",eval[,"beta"],"gamma",eval2[,"gamma"])),
                   nulling(dnorms_calc("delta",eval[,"delta"],"gamma",eval2[,"gamma"])), 
                   nulling(dnorms_calc("delta_minus",eval[,"delta_minus"],"gamma",eval2[,"gamma"])),
                   nulling(dnorms_calc("epsilon",eval[,"epsilon"],"gamma",eval2[,"gamma"])),
                   nulling(dnorms_calc("epsilon_minus",eval[,"epsilon_minus"],"gamma",eval2[,"gamma"])),
                   nulling(dnorms_calc("zeta",eval[,"zeta"],"gamma",eval2[,"gamma"])),
                   nulling(dnorms_calc("zeta_minus",eval[,"zeta_minus"],"gamma",eval2[,"gamma"])),
                   nulling(dnorms_calc("chi",eval[,"chi"],"gamma",eval2[,"gamma"])))
    
    delta_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"delta",eval2[,"delta"])),
                   nulling(dnorms_calc("beta",eval[,"beta"],"delta",eval2[,"delta"])),
                   nulling(dnorms_calc("gamma",eval[,"gamma"],"delta",eval2[,"delta"])),
                   nulling(dnorms_calc("epsilon",eval[,"epsilon"],"delta",eval2[,"delta"])),
                   nulling(dnorms_calc("zeta",eval[,"zeta"],"delta",eval2[,"delta"])),
                   nulling(dnorms_calc("chi",eval[,"chi"],"delta",eval2[,"delta"])))
    
    epsilon_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("alpha_plus",eval[,"alpha_plus"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("beta",eval[,"beta"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("beta_plus",eval[,"beta_plus"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("gamma",eval[,"gamma"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("gamma_plus",eval[,"gamma_plus"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("delta",eval[,"delta"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("delta_plus",eval[,"delta_plus"],"epsilon",eval2[,"epsilon"])),
                     nulling(dnorms_calc("zeta",eval[,"zeta"],"epsilon",eval2[,"epsilon"])))
    
    zeta_amb <- c(nulling(dnorms_calc("alpha_plus",eval[,"alpha_plus"],"zeta",eval2[,"zeta"])),
                  nulling(dnorms_calc("beta_plus",eval[,"beta_plus"],"zeta",eval2[,"zeta"])),
                  nulling(dnorms_calc("gamma_plus",eval[,"gamma_plus"],"zeta",eval2[,"zeta"])),
                  nulling(dnorms_calc("delta",eval[,"delta"],"zeta",eval2[,"zeta"])),
                  nulling(dnorms_calc("delta_plus",eval[,"delta_plus"],"zeta",eval2[,"zeta"])),
                  nulling(dnorms_calc("epsilon",eval[,"epsilon"],"zeta",eval2[,"zeta"])),
                  nulling(dnorms_calc("chi",eval[,"chi"],"zeta",eval2[,"zeta"])))
    
    chi_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"chi",eval2[,"chi"])),
                 nulling(dnorms_calc("beta",eval[,"beta"],"chi",eval2[,"chi"])),
                 nulling(dnorms_calc("gamma",eval[,"gamma"],"chi",eval2[,"chi"])),
                 nulling(dnorms_calc("delta",eval[,"delta"],"chi",eval2[,"chi"])),
                 nulling(dnorms_calc("epsilon",eval[,"epsilon"],"chi",eval2[,"chi"])),
                 nulling(dnorms_calc("zeta",eval[,"zeta"],"chi",eval2[,"chi"])))
    
    percent$alpha <- sum(alpha_amb) / length(alpha_amb)
    percent$beta <- sum(beta_amb) / length(beta_amb)
    percent$gamma <- sum(gamma_amb) / length(gamma_amb)
    percent$delta <- sum(delta_amb) / length(delta_amb)
    percent$epsilon <- sum(epsilon_amb) / length(epsilon_amb)
    percent$zeta <- sum(zeta_amb) / length(zeta_amb)
    percent$chi <- sum(chi_amb) / length(chi_amb)
    
    #Find how many 0s are found. 0s mean where original data angle is a NA
    ind_zero <- length(which(percent %in% 0))
    
    percent$total <- sum(c(percent$alpha, percent$beta, percent$gamma, percent$delta, percent$epsilon, percent$zeta, percent$chi)) / (7 - ind_zero)
    
    return(percent)
}



#Function to compute the dnorms
loopout <- function(ind = NULL, x = NULL,  eval2 = NULL){
    calcs <- c()
    for(i in 1:ind){
        calcs[i] <-  (dnorm(x = eval2, mean =  x$mean[i], sd = x$sd[i]) * x$prop[i]) / x$max_dens[i]
        #print(calcs)
        #print(calcs[i])
    }
    #print(calcs)
    calcs <- sum(calcs)
    return(calcs)
}

#Function to check over angle1 and angle 2 as well as the angle value
dnorms_calc <- function(angle = NULL, eval = NULL, angle2 = NULL, eval2 = NULL){
    if(is.na(eval) || is.na(eval2)){
        dnorming <- 0
    }else{
        ind_1 <- which(dict$total %in% angle)
        #let <- which(names(dict) %in% angle)
        ind_2 <- which(names(bd[[ind_1]][eval][[1]]) %in% angle2)
        
        #Extract data
        ind <- length(bd[[ind_1]][eval][[1]][[ind_2]]$prop)
        x <- bd[[ind_1]][eval][[1]][[ind_2]]
        dnorming <- loopout(ind = ind, x = x, eval2 = eval2 )
    }
    return(dnorming)
}


#Function to print 0 when we have Inf or NaN or NA
nulling <- function(x){
    if(length(x) == 0){
        x <- 0
    }else if(length(x) == 1){
        if (is.nan(x)){
            x <- 0
        }else if (is.na(x)){
            x <- 0
        }else if (x == Inf){
            x <- 0
        }
    }
    return(x)
}

#Normalization function
normalize <- function(data = NULL){
    max_v <- max(na.omit(data))
    min_v <- min(na.omit(data))
    
    data_norm <- c()
    for(i in 1:length(data)){
        data_norm[i] <- (data[i] - min_v) / (max_v - min_v)
    }
    
    return(data_norm)
}

#Create a new function to compute all dnorms for a specific file

neotest <- function(angle1 = NULL, eval = NULL, angle2=NULL, eval2 = NULL, bayesian_total = NULL){
    storing_values <- c()
    for(i in 1:100){
        bd <- bayesian_total[[i]]
        calc <- nulling(dnorms_calc(angle1, eval,angle2, eval2))
        storing_values[i] <- calc
    }
    
    #print(storing_values)
    store <- mean(storing_values)
    return(store)
}

## Replace pdb column o or b by given scores
replace_pdb <-
    function(score_vec, pdb, col="b") {
        nt_vec <- paste(pdb$atom$chain, pdb$atom$resno, pdb$atom$insert, sep="_")
        ntID_equivalences <- unique(nt_vec)
        
        ntID_vec <- unlist(lapply(1:nrow(pdb$atom), 
                                  function(i) which(nt_vec[i] == ntID_equivalences)))
        
        replacement <- score_vec[ntID_vec]
        
        pdb$atom[, col] <- replacement
        
        if (any(is.na(pdb$atom[, col]))) {
            ind <- which(is.na(pdb$atom[, col]))
            pdb$atom[ind, col] <- 0
        }
        
        return(pdb)
    }

