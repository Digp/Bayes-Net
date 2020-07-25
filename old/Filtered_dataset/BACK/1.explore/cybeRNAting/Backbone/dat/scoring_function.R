## Date: 2018-Feb-26 
#' Score function according with torsional angles data
#'
#' Function that takes a data.frame with the torsionals for every nucleotide
#' and computes score based on the reference data frame with the fitted 
#' gaussian models.
#'
#' @param torsionals A data.frame with the torsionals for all nucleotides,
#'     as obtained from measureNuc().
#' @param referencedata A data.frame with gaussians fitted to non-redundant
#'     RNA data.
#'
#' @return A data.frame with the scores per torsional and nucleotide
#'
#' @author Diego Gallego
#'
scoring_function <-
function(torsionals, referencedata) {
    scores <- lapply(torsionals$ntID,
                     .score_nt,
                     torsionals,
                     referencedata)
    out <- as.data.frame(matrix(unlist(scores), ncol=12, byrow=T))
    out <- cbind(torsionals$ntID, out)
    names(out) <- c("ntID", dihedrals)
    return(out)
}
.score_nt <-
function(ntID, torsionals, referencedata) {
    resid  <- torsionals[torsionals$ntID == ntID, "resid"]
    if (!(resid %in% c("A", "G", "C", "U"))) {
        return(rep(NA, 12))
    } else {
        .torsionals    <- torsionals[torsionals$ntID == ntID, ]
        .referencedata <- referencedata[referencedata$base == resid, ]
        scores <- lapply(dihedrals,
                        .score_tor,
                        .torsionals=.torsionals,
                        .referencedata=.referencedata)
        return(round(unlist(scores),3))
    }
}
.score_tor <-
function(tor, .torsionals, .referencedata) {
    ## get numeric value for given nt and angle
    angle <- .torsionals[, tor]

    if (is.na(angle)) {
        return(NA)

    } else {
        max_dens <- .referencedata[.referencedata$angle == tor, "max_dens"]
        mean     <- .referencedata[.referencedata$angle == tor, "mean"]
        sd       <- .referencedata[.referencedata$angle == tor, "sd"]
        prop     <- .referencedata[.referencedata$angle == tor, "weight"]

        s <- sum(dnorm(angle, mean, sd) * prop / max_dens)

        return(s)
    }
}
dihedrals <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi",
               "nu0", "nu1", "nu2", "nu3", "nu4")

##############################################################################

## Date: 2018-Feb-28
#' Take scores for each dihedral and nucleotide and compute backbone/sugar/chi
#' scores per nucleotide
#'
#' @param scores A data.frame as obtained from scoring_function()
#'
#' @return A vector with the averaged scores per nucleotide
#'
#' @author Diego Gallego
#'

score_backbone <-
function(scores) {
    return(round(
            rowMeans(
            scores[, c("alpha", "beta", "gamma", "delta", "epsilon", "zeta")], 
            na.rm=T), 2))
}

score_sugar <-
function(scores) {
    return(round(
            rowMeans(
            scores[, c("nu0", "nu1", "nu2", "nu3", "nu4")], 
            na.rm=T), 2))
}

score_chi <-
function(scores) {
    return(scores[, "chi"])
}

## Same idea, returns global score for the structure
score_str <-
function(scores) {
    return(round(mean(rowMeans(scores[, dihedrals], na.rm=T), na.rm=T), 3))
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

##############################################################################

## Date: 2018-Feb-28
#' Take scores for each dihedral and return outliers
#'
#' @param scores A data.frame as obtained from scoring_function()
#' @param threshold A numeric to select outliers
#' @param ntinfo A data.frame with information that relates the ntID with the 
#' residue identifier, chain and base.
#'
#' @return A data.frame
#'
#' @author Diego Gallego
#'

find_outliers <- 
function(scores, threshold, ntinfo) {
    outliers_ind <- which(scores[,dihedrals] < threshold, arr.ind=T)
    if (nrow(outliers_ind) == 0) {
        warning("No outliers found")
        return(NULL)
    }

    outliers_ind <- outliers_ind[order(outliers_ind[,1]),]
    if (any(ntinfo$insert != "?")) {
        ind <- which(ntinfo$insert != "?")
        ntinfo$resno[ind] <- paste(ntinfo$resno[ind], 
                                   ntinfo$insert[ind], sep="")
    }

    out <- data.frame(ntID=ntinfo$ntID[outliers_ind[,1]],
                      resid=ntinfo$resid[outliers_ind[,1]],
                      resno=ntinfo$resno[outliers_ind[,1]],
                      chain=ntinfo$chain[outliers_ind[,1]],
                      angle=dihedrals[outliers_ind[,2]], 
                      value=NA, score=NA, stringsAsFactors=F)

    for (i in 1:nrow(out)) {
        ntID <- out$ntID[i]
        out$value[i] <- ntinfo[ntinfo$ntID == ntID, out$angle[i]]
        out$score[i] <- scores[scores$ntID == ntID, out$angle[i]]
    }

    return(out[, -1])
}

