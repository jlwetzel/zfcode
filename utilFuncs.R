# Utility functions for manipulating b1hData, computing 
# log-probs, logOdds, regression models, etc ...

library(ggplot2)

# Some constants

caContacts <- list(c("n1", "a6"),
                   c("n2", "a3"),
                   c("n3", "a0"),
                   c("n3", "a2"))
bases <- c("A", "C", "G", "T")
residues <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
              "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# Some Functions

getFrame <- function(path) {
  # Convert the all.csv file at the specified path 
  # to a data frame and return it
  read.csv(paste0(path, "all.csv"))
}

getLogRatios <- function(dframe, countType = "actual",
                         scoreType = "relImpr",
                         contacts = caContacts) {
  # Compute the log-ratio scoring matrices for 
  # nuc-amino pairs in this dframe for the contacts 
  # specified and return them in a list.
  #
  # Two possible scoring ratios at present
  # 1. Relative improvement of binding probability:
  #  scoreType = "relImpr" --> lg(P(bind|N,R)/P(bind|N))
  # 2. Pairwise independence rating:
  #  scoreType = "indep"  --> lg( P(N,R|bind)/(P(N|bind)*P(R|bind)) )
  #
  # dframe - a dataframe with a column labeled "counts"
  #          and the remaining columns labeled as positions
  #          (e.g. n1, n2, ..., a1, a2, a3, ... ).
  # countType - decide whether to use raw counts ("normal") or 
  #             base 2 logs of counts ("log") when computing
  #             the empirical probs
  # scoreType - see above
  #
  # contacts - a list of 2d string vectors, one for each contact
  #            (e.g. see the caContacts constant above)

  if (scoreType == "relImpr")
    priorRatio <- 20  #P(N)/(P(N)*P(R))
  
  logOdds <- list()
  letterPairs <- as.vector(t(outer(bases, residues, 
                                   paste, sep = "")))
  
  # Create a log-ratio matrix for each contact
  for (con in contacts) {
    bpos <- con[1]
    apos <- con[2]
    contactLogOdds <- matrix(nrow = 4, ncol = 20)
    rownames(contactLogOdds) <- bases
    colnames(contactLogOdds) <- residues
    print(paste("contact:", bpos, apos))
    
    # Get the univariate and bivariate distributions
    # for this contact given binding
    bprobs <- getProbs(dframe, bpos, bases, countType)
    pairProbs <- getProbs(dframe, con, letterPairs, countType,
                          pairs = TRUE)
    if (scoreType == "indep")
      aprobs <- getProbs(dframe, apos, residues, countType)
    
    # Compute the log-ratios score for each pair at this contact
    for (base in bases)
      for (res in residues) {
        pair = paste0(base,res)

        if (scoreType == "relImpr")
        contactLogOdds[base,res] <- log2(pairProbs[pair]) -
                                    log2(bprobs[base]) + 
                                    log2(priorRatio)
        else if (scoreType == "indep")
        contactLogOdds[base,res] <- log2(pairProbs[pair]) -
                                    log2(bprobs[base] * aprobs[res]) 
      }
    logOdds[[paste0(con[1], con[2])]] <- contactLogOdds
  }
  logOdds
}

getProbs <- function(dframe, pos, letters, countType = "normal",
                     pairs = FALSE) {
  # Compute empirical probs from counts and +1 smoothing
  # Currently can compute empricial probabilities for any 
  # particular base position, residue position, or 
  # contact (i.e base and residue position pairing).
  #
  # dframe - a dataframe with a column labeled "counts"
  #          and the remaining columns labeled as positions
  #          (e.g. n1, n2, ..., a1, a2, a3, ... ).
  # pos - the position (or contact-pair) over which probs are to 
  #       be computed.  Should be a string giving the name of 
  #       the column in dframe.  In the case of a contact-pair, then
  #       pos should be a 2d vector of strings (e.g c("n3", "a0"))
  # letters - the complete alphabet over which the prob distribution
  #           is to be computed (as a string vector).  In the case 
  #           of a contact this should be the cross product of the 
  #           nuc and amino alphabets (i.e alphabet of all NR pairs)
  # countType - decide whether to use raw counts ("normal") or 
  #             base 2 logs of counts ("log") when computing
  #             the empirical probs
  # pairs - TRUE for contact-pairs, FALSE otherwise
  
  probs <- vector(mode = "numeric", length = length(letters))
  names(probs) <- letters
  
  if (countType == "actual")
    totCount <- sum(dframe$count)
  else
    totCount <- sum(floor(log2(dframe$count + 1)))

  # Perfrom the +1 smoothing
  for (letter in letters) {
    if (pairs) {
      base <- substr(letter, 1, 1)
      res <- substr(letter, 2, 2)
      count <- sum(dframe$count[dframe[pos[1]] == base
                                & dframe[pos[2]] == res])
      }
    else
      count <- sum(dframe$count[dframe[pos] == letter])
    if (count == 0)
      totCount <- totCount + 1
  }
  
  # Get the probs
  for (letter in letters) {
    if (pairs) {
      base <- substr(letter, 1, 1)
      res <- substr(letter, 2, 2)

      if (countType == "actual")
        count <- sum(dframe$count[dframe[pos[1]] == base 
                     & dframe[pos[2]] == res])
      else if (countType == "log")
        count <- sum(floor(log2(dframe$count[dframe[pos[1]] == base 
                     & dframe[pos[2]] == res] + 1)))
    }
    else {
      
      if (countType == "actual")
        count <- sum(dframe$count[dframe[pos] == letter])
      else if (countType == "log")
        count <- sum(floor(log2(dframe$count[dframe[pos] == 
                     letter] + 1)))
    }
    probs[letter] <- count/totCount
  }
  probs
}

getInfoContr(dframe, countType = "actual", contacts = caContacts) {
  # Returns a list of 2d matrices matrix of relative information 
  # contribution (one matrix per contact-pair) of eahc nuc-amino 
  # pairing in the joint distribution. 
}

scores2Heatmap <- function(matList){
  # Takes a list of 2D matrices of nuc-amino scores and 
  # returns a ggplot2 object that is a faceted heatmap 
  # with one panel for each matrix.

  # Convert the list of matrices to a dataframe 
  ratioFrame <- matList2frame(matList)

  p <- ggplot(ratioFrame, aes(x = rname, y = cname, 
              fill = val)) + 
       geom_tile() + facet_wrap(~key, nrow = 1) + 
       scale_fill_gradient2(name = "LogRatio",
                            low = "MediumOrchid3",
                            high = "ForestGreen",
                            mid = "grey97")  +
       xlab("Base") + 
       ylab("Residue") + 
       theme_bw()
  p
}

matList2frame <- function(matList) {
  # Takes a list of 2d matrices with rownames and 
  # columnnames (all same dimensions) and converts 
  # them to a corresponding dataframe
  # (i.e One observation for each entry of each matrix.
  #  The columns of the dataframe are the rowname, 
  #  colname, list-key, and actual value in the 
  #  matrix for the observation.)

  numObs <- (length(matList) * nrow(matList[[1]]) *
             ncol(matList[[1]]))
  print(numObs)
  vkey <- vector(mode = "character", length = numObs)
  vrname <- vector(mode = "character", length = numObs)
  vcname <- vector(mode = "character", length = numObs)
  vval <- vector(mode = "numeric", length = numObs)

  i <- 1
  for (key in names(matList)) 
    for (rname in rownames(matList[[key]]))
      for (cname in colnames(matList[[key]])) {
        vkey[i] <- key
        vrname[i] <- rname
        vcname[i] <- cname
        vval[i] <- matList[[key]][rname, cname]
        i = i + 1
      }
  data.frame(key = vkey, rname = vrname, cname = vcname,
             val = vval)
}