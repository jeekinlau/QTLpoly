
#' QTL allele effect estimation
#'
#' Computes allele specific and allele combination (within-parent) heritable effects from multiple QTL models.
#'
#' @param ploidy a numeric value of ploidy level of the cross (currently, only 2, 4 or 6).
#'
#' @param fitted a fitted multiple QTL model of class \code{qtlpoly.fitted}.
#'
#' @param x an object of class \code{qtlpoly.effects} to be plotted.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'fitted'} will be included.
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if \code{FALSE}, no output is produced.
#'
#' @param p1 a character string with the first parent name, e.g. \code{"P1"} (default).
#'
#' @param p2 a character string with the second parent name, e.g. \code{"P2"} (default).
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{qtlpoly.effects} which is a list of \code{results} for each containing the following components:
#'
#'     \item{pheno.col}{a phenotype column number.}
#'     \item{y.hat}{a vector with the predicted values.}
#'
#' @return A \pkg{ggplot2} barplot with parental allele and allele combination effects.
#'
#' @seealso \code{\link[qtlpoly]{read_data}}, \code{\link[qtlpoly]{remim}}, \code{\link[qtlpoly]{fit_model}}
#'
#' @examples
#'   \donttest{
#'   # Estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   library(qtlpoly)
#'   genoprob4x = lapply(maps4x[c(5)], calc_genoprob)
#'   data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1)
#'
#'   # Search for QTL
#'   remim.mod = remim(data = data, pheno.col = 1, w.size = 15, sig.fwd = 0.0011493379,
#' sig.bwd = 0.0002284465, d.sint = 1.5, n.clusters = 1)
#'
#'   # Fit model
#'   fitted.mod = fit_model(data, model=remim.mod, probs="joint", polygenes="none")
#'
#'   # Estimate effects
#'   est.effects = qtl_effects(ploidy = 4, fitted = fitted.mod, pheno.col = 1)
#'
#'   # Plot results
#'   plot(est.effects)
#'   }
#'   
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}, with modifications by Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \doi{10.1534/genetics.120.303080}.
#'     
#'     Kempthorne O (1955) The correlation between relatives in a simple autotetraploid population, \emph{Genetics} 40: 168-174.
#'
#' @export qtl_effects

qtl_effects <- function(ploidy = 6, fitted, pheno.col = NULL, verbose = TRUE) {
  
  if(is.null(pheno.col)) pheno.col <- fitted$pheno.col
  results <- vector("list", length(pheno.col))
  names(results) <- names(fitted$results)[which(fitted$pheno.col %in% pheno.col)]
  
  for(p in 1:length(results)) {
    
    if(!is.null(fitted$results[[names(results)[p]]]$qtls)) {
      
      nqtl <- dim(fitted$results[[names(results)[p]]]$qtls)[1]
      if(nqtl > 1) nqtl <- nqtl - 1
      effects <- vector("list", nqtl)
      qtl.mrk <- unlist(fitted$results[[names(results)[p]]]$qtls[c(1:nqtl),"Nmrk"])
      qtl.lg <- unlist(fitted$results[[names(results)[p]]]$qtls[c(1:nqtl),"LG"])
      qtl.pos <- unlist(fitted$results[[names(results)[p]]]$qtls[c(1:nqtl),"Pos"])
      
      if(verbose) {
        if(length(qtl.mrk) == 1) cat("There is ", length(qtl.mrk), " QTL in the model for trait ", pheno.col[p], " ", sQuote(names(results)[p]), ". Computing effects for QTL ", sep="")
        if(length(qtl.mrk) >= 2) cat("There are ", length(qtl.mrk), " QTL in the model for trait ", pheno.col[p], " ", sQuote(names(results)[p]), ". Computing effects for QTL ", sep="")
      }

      ## Creating temporary flag to indicate whether genotype probabilities come from MAPpoly or mappoly2
      if (length(fitted$results[[names(results)[p]]]$fitted$alleles) == ncol(combn(ploidy, ploidy/2))^2){
        genoprob_flag = TRUE
      } else {
        genoprob_flag = FALSE
      }

      if (!genoprob_flag) {
        for(q in 1:nqtl) {

          if(verbose) {
            if(q < nqtl) cat("...", qtl.mrk[q], "")
            if(q == nqtl) cat(paste0("... ", qtl.mrk[q]))
          }

          raw.blups <- fitted$results[[names(results)[p]]]$fitted$U[[q]]
          blups <- as.numeric(raw.blups)
          alleles <- as.character(fitted$results[[names(results)[p]]]$fitted$alleles)
          labels <- resolve_effect_labels(
            blups = raw.blups,
            alleles = alleles,
            fallback = fitted$homolog.names,
            context = paste0("trait ", pheno.col[p], ", QTL ", q)
          )
          names(blups) <- labels

          # For fit_model2 outputs, U already corresponds to homolog-level effects.
          effects[[q]] <- list(blups)
        }
      } else {
      
      if(ploidy == 6) {
        
        for(q in 1:nqtl) {
          
          if(verbose) {
            if(q < nqtl) cat("...", qtl.mrk[q], "")
            if(q == nqtl) cat(paste0("... ", qtl.mrk[q]))
          }
          
          blups <- fitted$results[[names(results)[p]]]$fitted$U[[q]]
          if (genoprob_flag == TRUE){
            alleles = matrix(unlist(strsplit(fitted$results[[names(results)[p]]]$fitted$alleles, '')), ncol=7, byrow=TRUE)[,-4]
            ## alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=7, byrow=TRUE)[,-4]
          } else {
            alleles = fitted$results[[names(results)[p]]]$fitted$alleles
          }
          
          A <- t(combn(letters[1:12],1))
          D <- t(combn(letters[1:12],2))
          T <- t(combn(letters[1:12],3))
          F <- t(combn(letters[1:12],4))
          G <- t(combn(letters[1:12],5))
          S <- t(combn(letters[1:12],6))
          
          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])
          t <- vector("list", dim(T)[1])
          f <- vector("list", dim(F)[1])
          g <- vector("list", dim(G)[1])
          s <- vector("list", dim(S)[1])

          if (genoprob_flag == TRUE){
            for(i in 1:dim(A)[1]) {
              a[[i]] <- which(alleles == as.character(A[i,1]), arr.ind = TRUE)[,1]
              a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
            }
            for(i in 1:dim(D)[1]) {
              d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
              d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
            }
            for(i in 1:dim(T)[1]) {
              t[[i]] <- which(apply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), 1, sum) == 3)
              t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
            }
            for(i in 1:dim(F)[1]) {
              f[[i]] <- which(apply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), 1, sum) == 4)
              f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
            }
            for(i in 1:dim(G)[1]) {
              g[[i]] <- which(apply(alleles == as.character(G[i,1]) | alleles == as.character(G[i,2]) | alleles == as.character(G[i,3]) | alleles == as.character(G[i,4]) | alleles == as.character(G[i,5]), 1, sum) == 5)
              g[[i]] <- mean(blups[Reduce(intersect, list(g[[i]]))])
            }
            for(i in 1:dim(S)[1]) {
              s[[i]] <- which(apply(alleles == as.character(S[i,1]) | alleles == as.character(S[i,2]) | alleles == as.character(S[i,3]) | alleles == as.character(S[i,4]) | alleles == as.character(S[i,5]) | alleles == as.character(S[i,6]), 1, sum) == 6)
              s[[i]] <- mean(blups[Reduce(intersect, list(s[[i]]))])
            }
          } else {
            for(i in 1:dim(A)[1]) {
              a[[i]] <- which(alleles == as.character(A[i,1]), arr.ind = TRUE)[1]
              a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
            }
            for(i in 1:dim(D)[1]) {
              d[[i]] <- which(lapply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), sum) == 1)
              d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
            }
            for(i in 1:dim(T)[1]) {
              t[[i]] <- which(lapply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), sum) == 1)
              t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
            }
            for(i in 1:dim(F)[1]) {
              f[[i]] <- which(lapply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), sum) == 1)
              f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
            }
            for(i in 1:dim(G)[1]) {
              g[[i]] <- which(lapply(alleles == as.character(G[i,1]) | alleles == as.character(G[i,2]) | alleles == as.character(G[i,3]) | alleles == as.character(G[i,4]) | alleles == as.character(G[i,5]), sum) == 1)
              g[[i]] <- mean(blups[Reduce(intersect, list(g[[i]]))])
            }
            for(i in 1:dim(S)[1]) {
              s[[i]] <- which(lapply(alleles == as.character(S[i,1]) | alleles == as.character(S[i,2]) | alleles == as.character(S[i,3]) | alleles == as.character(S[i,4]) | alleles == as.character(S[i,5]) | alleles == as.character(S[i,6]), sum) == 1)
              s[[i]] <- mean(blups[Reduce(intersect, list(s[[i]]))])
            }
          }
          names(a) <- as.character(A)
          names(d) <- apply(D, 1, paste, collapse="")
          names(t) <- apply(T, 1, paste, collapse="")
          names(f) <- apply(F, 1, paste, collapse="")
          names(g) <- apply(G, 1, paste, collapse="")
          names(s) <- apply(S, 1, paste, collapse="")
          
          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          t <- t[!is.nan(unlist(t))]
          f <- f[!is.nan(unlist(f))]
          g <- g[!is.nan(unlist(g))]
          s <- s[!is.nan(unlist(s))]
          
          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(t)) {
            t[[i]] <- t[[i]] -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(f)) {
            f[[i]] <- f[[i]] -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(g)) {
            g[[i]] <- g[[i]] -
              sum(unlist(f[which(lapply(lapply(strsplit(names(f), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 4) == TRUE)])) -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(g), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(s)) {
            s[[i]] <- s[[i]] -
              sum(unlist(g[which(lapply(lapply(strsplit(names(g), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 5) == TRUE)])) -
              sum(unlist(f[which(lapply(lapply(strsplit(names(f), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 4) == TRUE)])) -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(s), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          effects[[q]] <- list(unlist(a), unlist(d), unlist(t), unlist(f), unlist(g), unlist(s))
          
        }
      }
      
      if(ploidy == 4) {
        
        for(q in 1:nqtl) {
          
          if(verbose) {
            if(q < nqtl) cat("...", qtl.mrk[q], "")
            if(q == nqtl) cat(paste0("... ", qtl.mrk[q]))
          }
          
          blups <- fitted$results[[names(results)[p]]]$fitted$U[[q]]
          if (genoprob_flag == TRUE){
            alleles = matrix(unlist(strsplit(fitted$results[[names(results)[p]]]$fitted$alleles, '')), ncol=5, byrow=TRUE)[,-3]
            ## alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=5, byrow=TRUE)[,-3]
          } else {
            alleles = fitted$results[[names(results)[p]]]$fitted$alleles
          }
          
          A <- t(combn(letters[1:8],1))
          D <- t(combn(letters[1:8],2))
          T <- t(combn(letters[1:8],3))
          F <- t(combn(letters[1:8],4))
          
          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])
          t <- vector("list", dim(T)[1])
          f <- vector("list", dim(F)[1])

          if (genoprob_flag == TRUE){
            for(i in 1:dim(A)[1]) {
              a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[,1]
              a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
            }
            for(i in 1:dim(D)[1]) {
              d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
              d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
            }
            for(i in 1:dim(T)[1]) {
              t[[i]] <- which(apply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), 1, sum) == 3)
              t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
            }
            for(i in 1:dim(F)[1]) {
              f[[i]] <- which(apply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), 1, sum) == 4)
              f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
            }
          } else {
            for(i in 1:dim(A)[1]) {
              a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[1]
              a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
            } 
            for(i in 1:dim(D)[1]) {
              d[[i]] <- which(lapply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), sum) == 1)
              d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
            }
            for(i in 1:dim(T)[1]) {
              t[[i]] <- which(lapply(alleles == as.character(T[i,1]) | alleles == as.character(T[i,2]) | alleles == as.character(T[i,3]), sum) == 1)
              t[[i]] <- mean(blups[Reduce(intersect, list(t[[i]]))])
            }
            for(i in 1:dim(F)[1]) {
              f[[i]] <- which(lapply(alleles == as.character(F[i,1]) | alleles == as.character(F[i,2]) | alleles == as.character(F[i,3]) | alleles == as.character(F[i,4]), sum) == 1)
              f[[i]] <- mean(blups[Reduce(intersect, list(f[[i]]))])
            }
          }
          names(a) <- as.character(A)
          names(d) <- apply(D, 1, paste, collapse="")
          names(t) <- apply(T, 1, paste, collapse="")
          names(f) <- apply(F, 1, paste, collapse="")
          
          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          t <- t[!is.nan(unlist(t))]
          f <- f[!is.nan(unlist(f))]
          
          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(t)) {
            t[[i]] <- t[[i]] -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(t), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          for(i in 1:length(f)) {
            f[[i]] <- f[[i]] -
              sum(unlist(t[which(lapply(lapply(strsplit(names(t), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 3) == TRUE)])) -
              sum(unlist(d[which(lapply(lapply(strsplit(names(d), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 2) == TRUE)])) -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(f), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          effects[[q]] <- list(unlist(a), unlist(d), unlist(t), unlist(f))
          
        }
        
      }
      if(ploidy == 2) {
        
        for(q in 1:nqtl) {
          
          if(verbose) {
            if(q < nqtl) cat("...", qtl.mrk[q], "")
            if(q == nqtl) cat(paste0("... ", qtl.mrk[q]))
          }
          
          blups <- fitted$results[[names(results)[p]]]$fitted$U[[q]]

          ## Reassessing genoprob_flag for diploid case
          if (length(unlist(strsplit(fitted$results[[names(results)[p]]]$fitted$alleles, ''))) == 4) genoprob_flag = FALSE
          
          if (genoprob_flag == TRUE){
            alleles = matrix(unlist(strsplit(fitted$results[[names(results)[p]]]$fitted$alleles, '')), ncol=3, byrow=TRUE)[,-2]
            ## alleles <- matrix(unlist(strsplit(rownames(blups), '')), ncol=5, byrow=TRUE)[,-3]
          } else {
            alleles = fitted$results[[names(results)[p]]]$fitted$alleles
          }
          
          A <- t(combn(letters[1:4],1))
          D <- t(combn(letters[1:4],2))
          
          a <- vector("list", dim(A)[1])
          d <- vector("list", dim(D)[1])

          if (genoprob_flag == TRUE){
            for(i in 1:dim(A)[1]) {
              a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[,1]
              a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
            }
            for(i in 1:dim(D)[1]) {
              d[[i]] <- which(apply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), 1, sum) == 2)
              d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
            }
          } else {
            for(i in 1:dim(A)[1]) {
              a[[i]] <- which(alleles == as.character(A[i,]), arr.ind = TRUE)[1]
              a[[i]] <- mean(blups[Reduce(intersect, list(a[[i]]))])
            }
            for(i in 1:dim(D)[1]) {
              d[[i]] <- which(lapply(alleles == as.character(D[i,1]) | alleles == as.character(D[i,2]), sum) == 1)
              d[[i]] <- mean(blups[Reduce(intersect, list(d[[i]]))])
            }
          }
          names(a) <- as.character(A)
          names(d) <- apply(D, 1, paste, collapse="")
          
          a <- a[!is.nan(unlist(a))]
          d <- d[!is.nan(unlist(d))]
          
          for(i in 1:length(d)) {
            d[[i]] <- d[[i]] -
              sum(unlist(a[which(lapply(lapply(strsplit(names(a), split = ""), function(x) intersect(strsplit(names(d), split = "")[[i]], x)), function(x) length(x) == 1) == TRUE)]))
          }
          
          effects[[q]] <- list(unlist(a), unlist(d))
          
        }
        
      }
      }
      
      if(verbose) cat(". Done! \n\n", sep="")
      
    } else {
      
      if(verbose) cat("There are no QTL in the model for trait ", pheno.col[p], " ", sQuote(names(results)[p]), ". Skipping! \n\n", sep="")
      effects <- NULL
      qtl.lg <- NULL
      qtl.pos <- NULL
    }
    
    results[[p]] <- list(
      pheno.col=fitted$results[[names(results)[p]]]$pheno.col,
      qtl.lg=qtl.lg,
      qtl.pos=qtl.pos,
      effects=effects)
    
  }
  
  structure(list(fitted=deparse(substitute(fitted)),
                 ploidy=ploidy,
                 pheno.col=fitted$pheno.col,
          homolog.names=fitted$homolog.names,
          parent.names=fitted$parent.names,
                 results=results),
            class="qtlpoly.effects")
  
}

#' @rdname qtl_effects
#' @import ggplot2
#' @export

plot.qtlpoly.effects <- function(x, pheno.col = NULL, p1 = "P1", p2 = "P2", ...) {
  Alleles = Estimates = Parent = NULL
  parent.names <- x$parent.names
  homolog.names <- x$homolog.names
  if (identical(p1, "P1") && !is.null(parent.names) && length(parent.names) >= 1) {
    p1 <- parent.names[1]
  }
  if (identical(p2, "P2") && !is.null(parent.names) && length(parent.names) >= 2) {
    p2 <- parent.names[2]
  }
  if(is.null(pheno.col)) {
    pheno.col <- seq_along(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }

  infer_parent <- function(alleles, parent.names = NULL, homolog.names = NULL) {
    alleles <- as.character(alleles)
    parent.pattern <- "^(.+?)[_.-]h[0-9]+$"

    # Case 1: labels already encode parent names (e.g. Parent_h1).
    if (all(grepl(parent.pattern, alleles))) {
      return(sub(parent.pattern, "\\1", alleles))
    }

    # Case 2: labels match homolog names that encode parent names.
    if (!is.null(homolog.names) && length(homolog.names) > 0) {
      homolog.names <- as.character(homolog.names)
      if (!anyDuplicated(homolog.names) &&
          all(alleles %in% homolog.names) &&
          all(grepl(parent.pattern, homolog.names))) {
        idx <- match(alleles, homolog.names)
        return(sub(parent.pattern, "\\1", homolog.names[idx]))
      }
    }

    # Case 3: conservative fallback for biparental blocks where homolog labels
    # do not encode parent names (e.g. h1..h8 or a..h).
    if (!is.null(parent.names) && length(parent.names) > 0) {
      parent.names <- unique(as.character(parent.names[!is.na(parent.names) & parent.names != ""]))
      if (length(parent.names) > 2) {
        stop(
          "Could not infer parent names from allele labels in a multiparent setting. ",
          "Please use allele/homolog labels in the format '<Parent>_h<number>'."
        )
      }

      if (length(parent.names) != 2) {
        stop(
          "Could not infer parent names from allele labels. ",
          "Biparental fallback requires exactly two parent names."
        )
      }

      ref.homologs <- homolog.names
      if (is.null(ref.homologs) || length(ref.homologs) == 0) {
        ref.homologs <- alleles
      }
      ref.homologs <- as.character(ref.homologs)

      if (!anyDuplicated(ref.homologs) &&
          length(alleles) == length(ref.homologs) &&
          length(ref.homologs) %% length(parent.names) == 0) {
        idx <- match(alleles, ref.homologs)
        if (!any(is.na(idx))) {
          homologs.per.parent <- length(ref.homologs) / length(parent.names)
          parent.idx <- ((idx - 1) %/% homologs.per.parent) + 1
          if (all(parent.idx >= 1 & parent.idx <= length(parent.names))) {
            return(parent.names[parent.idx])
          }
        }
      }
    }

    stop(
      "Could not infer parent names from allele labels. ",
      "Expected '<Parent>_h<number>' labels, or biparental labels that match ordered homolog metadata."
    )
  }

  res = list()
  for(p in pheno.col) {
    nqtl <- length(x$results[[p]]$effects)
    qtl.lg <- x$results[[p]]$qtl.lg
    qtl.pos <- x$results[[p]]$qtl.pos
    if(nqtl > 0) {
      for(q in 1:nqtl) {
        add <- x$results[[p]]$effects[[q]][[1]]
        if (is.null(add)) next
        add <- as.numeric(add)
        add.names <- names(x$results[[p]]$effects[[q]][[1]])
        if (is.null(add.names) || length(add.names) != length(add) || any(is.na(add.names)) || any(add.names == "")) {
          stop("Missing or invalid allele labels for trait ", names(x$results)[p], ", QTL ", q, ".")
        }
        if (anyDuplicated(add.names)) {
          stop("Duplicated allele labels for trait ", names(x$results)[p], ", QTL ", q, ".")
        }

        data <- data.frame(
          Estimates = add,
          Alleles = add.names,
          stringsAsFactors = FALSE
        )
        data$Parent <- infer_parent(
          alleles = data$Alleles,
          parent.names = parent.names,
          homolog.names = homolog.names
        )
        if (!is.null(parent.names) && length(parent.names) > 0) {
          parent.levels <- unique(c(parent.names, unique(data$Parent)))
        } else {
          parent.levels <- unique(data$Parent)
        }
        data$Parent <- factor(data$Parent, levels = parent.levels)

        subtitle.text <- paste("QTL", q)
        if (!is.null(qtl.lg) && length(qtl.lg) >= q && !is.na(qtl.lg[q]) && qtl.lg[q] != "") {
          subtitle.text <- paste0("QTL ", q, ": LG ", qtl.lg[q])
        }
        if (!is.null(qtl.pos) && length(qtl.pos) >= q && !is.na(qtl.pos[q]) && qtl.pos[q] != "") {
          pos.value <- as.character(qtl.pos[q])
          pos.num <- suppressWarnings(as.numeric(pos.value))
          if (!is.na(pos.num)) {
            pos.value <- format(round(pos.num, 2), nsmall = 0, trim = TRUE, scientific = FALSE)
          }
          if (!grepl("cM$", pos.value, ignore.case = TRUE)) {
            pos.value <- paste0(pos.value, " cM")
          }
          subtitle.text <- paste0(subtitle.text, " at ", pos.value)
        }

        plot <- ggplot(data, aes(x = Alleles, y = Estimates, fill = Estimates)) +
          geom_bar(stat="identity") +
          scale_fill_gradient2(low = "red", high = "blue", guide = "none") +
          labs(title=names(x$results)[p], subtitle=paste0(subtitle.text, "\n")) +
          facet_wrap(. ~ Parent, scales="free_x", ncol = min(4, length(unique(data$Parent))), strip.position="bottom") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.x.bottom = element_text(hjust = 1, vjust = 0.5))
        res = c(res, plot)
        print(plot)
      }
    }
  }
  invisible(res)
}

