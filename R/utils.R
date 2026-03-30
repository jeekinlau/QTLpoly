#' @keywords internal
last <- function(x) { tail(x, n = 1) }

#' @keywords internal
h <- function(w) if( any( grepl( "Recycling array of length 1 in vector-array arithmetic is deprecated", w) ) ) invokeRestart( "muffleWarning" )

#' @keywords internal
sanitize_homolog_names <- function(homolog.names, n.levels, prefix = "h") {
	if (is.null(n.levels) || n.levels < 1) return(character(0))
	if (is.null(homolog.names) || length(homolog.names) != n.levels) {
		homolog.names <- paste0(prefix, seq_len(n.levels))
	}
	homolog.names <- as.character(homolog.names)
	bad <- is.na(homolog.names) | homolog.names == ""
	if (any(bad)) {
		homolog.names[bad] <- paste0(prefix, which(bad))
	}
	make.unique(homolog.names, sep = "_dup")
}

#' @keywords internal
infer_parent_names_from_homologs <- function(homolog.names) {
	if (is.null(homolog.names) || length(homolog.names) == 0) return(NULL)
	if (all(grepl("_h[0-9]+$", homolog.names))) {
		return(unique(sub("_h[0-9]+$", "", homolog.names)))
	}
	NULL
}

#' @keywords internal
normalize_homolog_parent_metadata <- function(Z, parent.names = NULL) {
	if (is.null(Z)) {
		return(list(Z = Z, homolog.names = NULL, parent.names = parent.names))
	}
	homolog.names <- sanitize_homolog_names(dimnames(Z)[[1]], nrow(Z))
	dimnames(Z)[[1]] <- homolog.names

	if (is.null(parent.names) || length(parent.names) == 0 || any(is.na(parent.names)) || any(parent.names == "")) {
		parent.names <- infer_parent_names_from_homologs(homolog.names)
	}
	if (!is.null(parent.names)) {
		parent.names <- unique(as.character(parent.names[!is.na(parent.names) & parent.names != ""]))
		if (length(parent.names) == 0) {
			parent.names <- NULL
		}
	}

	list(Z = Z, homolog.names = homolog.names, parent.names = parent.names)
}

#' @keywords internal
align_random_effect_names <- function(Z, K, effect.name = "random effect") {
	if (is.null(Z) || is.null(K)) {
		return(list(Z = Z, K = K))
	}
	if (!is.matrix(K)) {
		K <- as.matrix(K)
	}
	if (nrow(K) != ncol(K)) {
		stop("Variance-covariance matrix K for ", effect.name, " must be square.")
	}

	z.names <- colnames(Z)
	if (is.null(z.names) || length(z.names) != ncol(Z) || any(is.na(z.names)) || any(z.names == "") || anyDuplicated(z.names)) {
		k.seed <- rownames(K)
		if (is.null(k.seed) || length(k.seed) != ncol(Z) || any(is.na(k.seed)) || any(k.seed == "") || anyDuplicated(k.seed)) {
			z.names <- paste0("lvl", seq_len(ncol(Z)))
		} else {
			z.names <- as.character(k.seed)
		}
		colnames(Z) <- z.names
	} else {
		z.names <- as.character(z.names)
		colnames(Z) <- z.names
	}

	if (is.null(rownames(K)) && !is.null(colnames(K))) {
		rownames(K) <- colnames(K)
	}
	if (is.null(colnames(K)) && !is.null(rownames(K))) {
		colnames(K) <- rownames(K)
	}
	if (is.null(rownames(K)) || is.null(colnames(K))) {
		if (ncol(Z) != ncol(K)) {
			stop("Dimension mismatch between Z and K for ", effect.name, ".")
		}
		rownames(K) <- colnames(K) <- z.names
	}

	k.row <- as.character(rownames(K))
	k.col <- as.character(colnames(K))
	if (anyDuplicated(k.row) || anyDuplicated(k.col)) {
		stop("Duplicated level names in K for ", effect.name, ".")
	}
	if (!setequal(k.row, k.col)) {
		stop("Row and column names of K differ for ", effect.name, ".")
	}
	if (ncol(Z) != nrow(K)) {
		stop("Dimension mismatch between Z and K for ", effect.name, ".")
	}
	if (!setequal(z.names, k.row)) {
		missing.in.k <- setdiff(z.names, k.row)
		missing.in.z <- setdiff(k.row, z.names)
		stop(
			"Name mismatch between Z and K for ", effect.name,
			". Missing in K: ", paste(missing.in.k, collapse = ", "),
			". Missing in Z: ", paste(missing.in.z, collapse = ", "),
			"."
		)
	}

	K <- K[z.names, z.names, drop = FALSE]
	list(Z = Z, K = K)
}

#' @keywords internal
align_eta_names <- function(ETA, eta.label = "ETA") {
	if (is.null(ETA)) {
		return(ETA)
	}
	if (!is.list(ETA)) {
		stop(eta.label, " must be a list.")
	}

	for (i in seq_along(ETA)) {
		if (!is.list(ETA[[i]]) || !all(c("Z", "K") %in% names(ETA[[i]]))) {
			next
		}
		effect.name <- names(ETA)[i]
		if (is.null(effect.name) || effect.name == "") {
			effect.name <- paste0(eta.label, "[", i, "]")
		}
		aligned <- align_random_effect_names(Z = ETA[[i]]$Z, K = ETA[[i]]$K, effect.name = effect.name)
		ETA[[i]]$Z <- aligned$Z
		ETA[[i]]$K <- aligned$K
	}
	ETA
}

#' @keywords internal
resolve_effect_labels <- function(blups, alleles = NULL, fallback = NULL, context = "effect") {
	n <- length(blups)
	valid_labels <- function(x) {
		!is.null(x) && length(x) == n && !any(is.na(x)) && !any(x == "") && !anyDuplicated(x)
	}

	labels <- NULL
	if (valid_labels(names(blups))) {
		labels <- as.character(names(blups))
	}
	if (is.null(labels) && valid_labels(alleles)) {
		labels <- as.character(alleles)
	}
	if (is.null(labels) && valid_labels(fallback)) {
		labels <- as.character(fallback)
	}
	if (is.null(labels)) {
		stop("Could not infer unique homolog labels for ", context, ".")
	}

	if (valid_labels(alleles) && !identical(labels, as.character(alleles))) {
		stop("Homolog label order mismatch detected for ", context, ".")
	}
	labels
}
