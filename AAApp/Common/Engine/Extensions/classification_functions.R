#' Classification functions for Advanced_Amplicon_Analysis
#'
#' Each function accepts a named logical vector and returns one category.

classify_ectoine <- function(x) {
  ect <- as.logical(x[c("ectA", "ectB", "ectC")])
  ectD <- as.logical(x["ectD"])
  if (all(is.na(c(ect, ectD)))) {
    return("Unknown ectoine biosynthetic potential")
  }
  n <- sum(ect %in% TRUE)
  if (n == 3 && isTRUE(ectD)) {
    return("Ectoine and hydroxyectoine producer")
  }
  if (n == 3) {
    return("Ectoine producer")
  }
  if (n == 2) {
    return("High ectoine biosynthetic potential")
  }
  if (n == 1) {
    return("Partial ectoine biosynthetic potential")
  }
  "No ectoine biosynthetic potential"
}

classify_methanotroph <- function(x) {
  pmo <- as.logical(x[c("pmoA", "pmoB", "pmoC")])
  mmo <- as.logical(x[c("mmoX", "mmoY", "mmoZ", "mmoB", "mmoC", "mmoD")])
  if (all(is.na(c(pmo, mmo)))) {
    return("Unknown methanotrophic potential")
  }
  p <- sum(pmo %in% TRUE)
  m <- sum(mmo %in% TRUE)
  if (p == 3 && m == 6) {
    return("pMMO + sMMO methanotroph")
  }
  if (p == 3) {
    return("Complete pMMO methanotroph")
  }
  if (m == 6) {
    return("Complete sMMO methanotroph")
  }
  if (p >= 2) {
    return("Probable pMMO methanotroph")
  }
  if (m >= 3) {
    return("Probable sMMO methanotroph")
  }
  if (p >= 1 || m >= 1) {
    return("Possible methanotroph")
  }
  "Non-methanotroph"
}

classify_methanogenesis <- function(x) {
  mcr <- as.logical(x[c("mcrA", "mcrB", "mcrG")])
  if (all(is.na(mcr))) {
    return("Unknown methanogenic potential")
  }
  n <- sum(mcr %in% TRUE)
  if (n == 3) {
    return("Methanogenic potential")
  }
  if (n >= 1) {
    return("Partial methanogenic potential")
  }
  "Non-methanogen"
}

classify_nitrification <- function(x) {
  amo <- as.logical(x[c("amoA", "amoB", "amoC")])
  hao <- as.logical(x["hao"])
  nxr <- as.logical(x[c("nxrA", "nxrB")])
  if (all(is.na(c(amo, hao, nxr)))) {
    return("Unknown nitrification potential")
  }
  amo_complete <- all(amo %in% TRUE)
  nxr_complete <- all(nxr %in% TRUE)
  if (amo_complete && isTRUE(hao) && nxr_complete) {
    return("Complete nitrifier")
  }
  if (amo_complete && isTRUE(hao)) {
    return("Ammonia oxidizer")
  }
  if (nxr_complete) {
    return("Nitrite oxidizer")
  }
  if (any(c(amo, hao, nxr) %in% TRUE)) {
    return("Partial nitrification potential")
  }
  "Non-nitrifier"
}

classify_denitrification <- function(x) {
  genes <- as.logical(x[c("narG", "nirK", "nirS", "norB", "nosZ")])
  if (all(is.na(genes))) {
    return("Unknown denitrification potential")
  }
  n <- sum(genes %in% TRUE)
  if (n == length(genes)) {
    return("Complete denitrifier")
  }
  if (n >= 3) {
    return("Partial denitrifier")
  }
  if (n > 0) {
    return("Limited denitrification potential")
  }
  "Non-denitrifier"
}

classify_sulfate_reducer <- function(x) {
  genes <- as.logical(x[c("sat", "aprA", "aprB", "dsrA", "dsrB")])
  if (all(is.na(genes))) {
    return("Unknown sulfate-reduction potential")
  }
  if (all(genes %in% TRUE)) {
    return("Sulfate reducer")
  }
  if (any(genes[4:5] %in% TRUE)) {
    return("Potential sulfate reducer")
  }
  "Non-sulfate reducer"
}

classify_nitrogen_fixation <- function(x) {
  genes <- as.logical(x[c("nifH", "nifD", "nifK")])
  if (all(is.na(genes))) {
    return("Unknown nitrogen-fixation potential")
  }
  if (all(genes %in% TRUE)) {
    return("Nitrogen fixer")
  }
  if (any(genes %in% TRUE)) {
    return("Potential nitrogen fixer")
  }
  "Non-nitrogen fixer"
}

classify_anammox <- function(x) {
  genes <- as.logical(x[c("hzsA", "hzsB", "hzsC", "hdh")])
  if (all(is.na(genes))) {
    return("Unknown anammox potential")
  }
  if (all(genes %in% TRUE)) {
    return("Anammox bacterium")
  }
  if (any(genes %in% TRUE)) {
    return("Potential anammox bacterium")
  }
  "Non-anammox bacterium"
}

classify_acetogenesis <- function(x) {
  genes <- as.logical(x[c("acsA", "acsB", "fhs", "fdh")])
  if (all(is.na(genes))) {
    return("Unknown acetogenic potential")
  }
  if (all(genes %in% TRUE)) {
    return("Acetogen")
  }
  if (any(genes %in% TRUE)) {
    return("Potential acetogen")
  }
  "Non-acetogen"
}


#' Normalise named gene-presence values used by C1 classifiers
#'
#' Accepts logical, numeric and common textual presence/absence encodings.
aaa_normalise_gene_presence <- function(x) {
  if (is.null(names(x))) stop("x must be a named vector.", call. = FALSE)
  out <- rep(NA, length(x))
  names(out) <- names(x)
  if (is.logical(x)) {
    return(x)
  }
  if (is.numeric(x)) {
    out[!is.na(x)] <- x[!is.na(x)] > 0
    return(as.logical(out))
  }
  value <- trimws(tolower(as.character(x)))
  out[value %in% c("true", "t", "yes", "y", "present", "detected", "1")] <- TRUE
  out[value %in% c("false", "f", "no", "n", "absent", "not detected", "0")] <- FALSE
  as.logical(out)
}

#' Evaluate methanogenesis and reductive Wood-Ljungdahl evidence
#'
#' Internal structured evaluator shared by the combined and route-specific
#' classifiers. Methanogenesis has interpretative priority over
#' homoacetogenesis: a reductive WLP-like gene set in a genome with a supported
#' Mcr terminal complex is not labelled homoacetogenic.
aaa_evaluate_methanogenesis_wlp <- function(x) {
  x <- aaa_normalise_gene_presence(x)

  get_values <- function(genes) {
    values <- stats::setNames(rep(NA, length(genes)), genes)
    shared <- intersect(names(x), genes)
    values[shared] <- x[shared]
    values
  }
  present <- function(genes) {
    values <- get_values(genes)
    stats::setNames(values %in% TRUE, genes)
  }
  count_present <- function(genes) sum(present(genes))
  any_present <- function(genes) any(present(genes))
  all_present <- function(genes) all(present(genes))

  terminal <- c("mcrA", "mcrB", "mcrG")
  co2_core <- c("ftr", "mch", "mtd", "mer")
  mtr <- paste0("mtr", LETTERS[1:8])
  acetate_activation <- c("ackA", "pta", "acs")
  cdh <- paste0("cdh", LETTERS[1:5])
  methanol <- c("mtaA", "mtaB", "mtaC")
  methylamine_modules <- list(
    monomethylamine = c("mtmB", "mtmC"),
    dimethylamine = c("mtbB", "mtbC"),
    trimethylamine = c("mttB", "mttC")
  )
  methylamine_all <- unique(c("mtbA", unlist(methylamine_modules)))
  wlp_carbonyl <- c("acsA", "acsB")
  wlp_corrinoid <- c("acsC", "acsD", "acsE")
  wlp_methyl <- c("fhs", "folD", "metF")
  # 'fdhA' was declared by no function, so it was never searched; the fdh
  # dictionary entry already aliases "formate dehydrogenase subunit alpha".
  wlp_electron <- c("fdh")

  requested <- unique(c(
    terminal, co2_core, mtr, acetate_activation, cdh, methanol,
    methylamine_all, wlp_carbonyl, wlp_corrinoid, wlp_methyl, wlp_electron,
    paste0("fwd", LETTERS[1:4]), paste0("fmd", LETTERS[1:4])
  ))
  unknown <- all(is.na(get_values(requested)))

  terminal_complete <- all_present(terminal)
  terminal_supported <- isTRUE(present("mcrA")[[1]]) &&
    any_present(c("mcrB", "mcrG"))

  fwd_fmd_subunits <- vapply(
    LETTERS[1:4],
    function(subunit) any_present(c(paste0("fwd", subunit), paste0("fmd", subunit))),
    logical(1)
  )

  hydrogen_complete <- terminal_supported &&
    count_present(co2_core) >= 3L &&
    sum(fwd_fmd_subunits) >= 3L &&
    count_present(mtr) >= 6L

  acetate_activation_complete <- all_present(c("ackA", "pta")) ||
    isTRUE(present("acs")[[1]])
  acetoclastic_complete <- terminal_supported &&
    acetate_activation_complete && count_present(cdh) >= 4L

  methanol_entry_complete <- count_present(methanol) >= 2L
  methylamine_entry_complete <- any(vapply(
    methylamine_modules,
    function(module) all_present(module),
    logical(1)
  ))
  methylotrophic_complete <- terminal_supported &&
    (methanol_entry_complete || methylamine_entry_complete)

  wlp_complete <- all_present(wlp_carbonyl) &&
    count_present(wlp_corrinoid) >= 2L &&
    count_present(wlp_methyl) >= 2L &&
    any_present(wlp_electron)

  methanogenic_routes <- character()
  if (hydrogen_complete) methanogenic_routes <- c(methanogenic_routes, "Hydrogenotrophic methanogenesis")
  if (acetoclastic_complete) methanogenic_routes <- c(methanogenic_routes, "Acetoclastic methanogenesis")
  if (methylotrophic_complete) methanogenic_routes <- c(methanogenic_routes, "Methylotrophic methanogenesis")

  relevant_marker_detected <- any_present(requested)

  list(
    unknown = unknown,
    terminal_complete = terminal_complete,
    terminal_supported = terminal_supported,
    hydrogen_complete = hydrogen_complete,
    acetoclastic_complete = acetoclastic_complete,
    methylotrophic_complete = methylotrophic_complete,
    wlp_complete = wlp_complete,
    methanogenic_routes = methanogenic_routes,
    relevant_marker_detected = relevant_marker_detected
  )
}

#' Classify methanogenesis subtypes and non-methanogenic WLP potential
#'
#' Complete methanogenic routes are reported first. Homoacetogenesis is assigned
#' only when a complete reductive Wood-Ljungdahl profile is detected without a
#' supported methanogenic terminal complex. Partial profiles are grouped as
#' insufficient evidence rather than labelled probable.
classify_methanogenesis_subtype <- function(x) {
  ev <- aaa_evaluate_methanogenesis_wlp(x)
  if (ev$unknown) {
    return("Unknown methanogenesis and Wood-Ljungdahl profile")
  }
  if (length(ev$methanogenic_routes)) {
    return(paste0(paste(ev$methanogenic_routes, collapse = " + "), " potential"))
  }
  if (!ev$terminal_supported && ev$wlp_complete) {
    return("Reductive Wood-Ljungdahl pathway potential (possible homoacetogenesis)")
  }
  if (ev$relevant_marker_detected) {
    return("Insufficient evidence to classify metabolic potential")
  }
  "No detected methanogenesis or Wood-Ljungdahl pathway potential"
}

#' Classify hydrogenotrophic methanogenesis potential
classify_hydrogenotrophic_methanogenesis <- function(x) {
  ev <- aaa_evaluate_methanogenesis_wlp(x)
  if (ev$unknown) {
    return("Unknown hydrogenotrophic methanogenesis potential")
  }
  if (ev$hydrogen_complete) {
    return("Hydrogenotrophic methanogenesis potential")
  }
  if (ev$relevant_marker_detected) {
    return("Insufficient evidence to classify hydrogenotrophic methanogenesis potential")
  }
  "No detected hydrogenotrophic methanogenesis potential"
}

#' Classify acetoclastic methanogenesis potential
classify_acetoclastic_methanogenesis <- function(x) {
  ev <- aaa_evaluate_methanogenesis_wlp(x)
  if (ev$unknown) {
    return("Unknown acetoclastic methanogenesis potential")
  }
  if (ev$acetoclastic_complete) {
    return("Acetoclastic methanogenesis potential")
  }
  if (ev$relevant_marker_detected) {
    return("Insufficient evidence to classify acetoclastic methanogenesis potential")
  }
  "No detected acetoclastic methanogenesis potential"
}

#' Classify methylotrophic methanogenesis potential
#'
#' Requires a supported Mcr terminal complex plus a coherent methanol,
#' monomethylamine, dimethylamine or trimethylamine substrate-entry module.
classify_methylotrophic_methanogenesis <- function(x) {
  ev <- aaa_evaluate_methanogenesis_wlp(x)
  if (ev$unknown) {
    return("Unknown methylotrophic methanogenesis potential")
  }
  if (ev$methylotrophic_complete) {
    return("Methylotrophic methanogenesis potential")
  }
  if (ev$relevant_marker_detected) {
    return("Insufficient evidence to classify methylotrophic methanogenesis potential")
  }
  "No detected methylotrophic methanogenesis potential"
}

#' Score reductive Wood-Ljungdahl pathway evidence
#'
#' The score measures WLP gene-module completeness, not a demonstrated
#' homoacetogenic phenotype. Methanogenic genomes are explicitly excluded from
#' homoacetogenic interpretation while retaining the pathway-completeness score.
score_homoacetogenesis <- function(x) {
  required_names <- c(
    "acsA", "acsB", "acsC", "acsD", "acsE",
    "fhs", "folD", "metF", "fdh", "pta", "ackA",
    "mcrA", "mcrB", "mcrG"
  )
  values <- stats::setNames(rep(NA, length(required_names)), required_names)
  shared <- intersect(names(x), required_names)
  values[shared] <- aaa_normalise_gene_presence(x[shared])
  if (all(is.na(values))) {
    return(list(score = NA_real_, level = "Unknown", evaluated_fraction = 0))
  }

  present <- stats::setNames(values %in% TRUE, names(values))
  evaluated <- stats::setNames(!is.na(values), names(values))
  group_score <- function(genes, weight, minimum = 1L) {
    if (!any(evaluated[genes])) {
      return(0)
    }
    weight * min(1, sum(present[genes]) / minimum)
  }
  score <- group_score(c("acsA", "acsB"), 40, 2) +
    group_score(c("acsC", "acsD", "acsE"), 25, 2) +
    group_score(c("fhs", "folD", "metF"), 20, 2) +
    group_score(c("fdh"), 10, 1) +
    group_score(c("pta", "ackA"), 5, 2)
  score <- round(min(100, score), 1)
  methanogenic_exclusion <- isTRUE(present[["mcrA"]]) &&
    any(present[c("mcrB", "mcrG")])
  level <- if (methanogenic_exclusion) "Excluded: methanogenic pathway detected" else if (score >= 80) "High" else if (score > 0) "Insufficient" else "None"
  list(
    score = score,
    level = level,
    evaluated_fraction = round(sum(evaluated) / length(values), 3),
    methanogenic_exclusion = methanogenic_exclusion
  )
}

#' Classify homoacetogenic potential
#'
#' A complete reductive WLP profile is interpreted as possible
#' homoacetogenesis only in the absence of a supported Mcr terminal complex.
classify_homoacetogenesis <- function(x) {
  ev <- aaa_evaluate_methanogenesis_wlp(x)
  if (ev$unknown) {
    return("Unknown homoacetogenic potential")
  }
  if (ev$terminal_supported) {
    return("No homoacetogenesis assigned: methanogenic pathway detected")
  }
  if (ev$wlp_complete) {
    return("Homoacetogenesis potential")
  }
  if (ev$relevant_marker_detected) {
    return("Insufficient evidence to classify homoacetogenic potential")
  }
  "No detected homoacetogenic potential"
}

# =============================================================================
# Shared functional-abundance selector.
#
# Every curated function exposes a `pathway_selector` that decides which
# classifier outputs count as a positive functional call for the potential
# metabolomic pathways abundance analysis. This single predicate replaces the
# three ad-hoc idioms that had drifted apart (`^(No detected|Unknown)`,
# `^(Non-|Unknown)` plus a keyword, and hard-coded string vectors). It treats a
# call as positive unless it is a negative, unknown, insufficient or excluded
# statement, so partial-but-real evidence still counts while "Insufficient
# evidence…" and "No … assigned" no longer leak through as positives.
#
# Two functions intentionally keep a stricter, positive-keyword selector
# (Syngas_to_ethanol requires a resolved ethanol branch; PLA_engineering_chassis
# requires a "Putative" call), so they are deliberately NOT migrated here.
aaa_positive_pathway_call <- function(x) {
  !grepl("^(no |non-|insufficient|unknown|excluded)", x, ignore.case = TRUE)
}

# =============================================================================
# Generic pathway classifier used by curated biological functions.
# It keeps the existing genome/GFF search engine unchanged and standardises
# Complete / Probable / Partial / No detected / Unknown decisions.
# =============================================================================
aaa_make_pathway_classifier <- function(label,
                                        required_groups,
                                        supporting_groups = list(),
                                        complete_required_fraction = 1,
                                        probable_required_fraction = 0.5) {
  force(label)
  force(required_groups)
  force(supporting_groups)
  function(x) {
    original_names <- names(x)
    x <- as.logical(x)
    names(x) <- original_names
    requested <- unique(c(
      unlist(required_groups, use.names = FALSE),
      unlist(supporting_groups, use.names = FALSE)
    ))
    if (length(requested) == 0L || all(is.na(x[requested]))) {
      return(paste("Unknown", label, "potential"))
    }
    present <- function(genes) {
      vals <- x[genes]
      vals[is.na(vals)] <- FALSE
      sum(vals %in% TRUE)
    }
    group_fraction <- function(groups) {
      if (!length(groups)) {
        return(numeric())
      }
      vapply(groups, function(g) present(g) / max(1L, length(g)), numeric(1))
    }
    req <- group_fraction(required_groups)
    sup <- group_fraction(supporting_groups)
    complete <- length(req) > 0L && all(req >= complete_required_fraction)
    probable <- length(req) > 0L && all(req >= probable_required_fraction) &&
      (length(sup) == 0L || any(sup > 0))
    partial <- any(c(req, sup) > 0)
    if (complete) {
      return(paste("Complete", label, "potential"))
    }
    if (probable) {
      return(paste("Probable", label, "potential"))
    }
    if (partial) {
      return(paste("Partial", label, "evidence"))
    }
    paste("No detected", label, "potential")
  }
}

# =============================================================================
# Additional non-duplicated functional classifiers
# Pigments, compatible-solute precursors and C1/syngas valorisation.
# =============================================================================

aaa_marker_vector <- function(x, genes) {
  values <- stats::setNames(rep(NA, length(genes)), genes)
  shared <- intersect(names(x), genes)
  values[shared] <- as.logical(x[shared])
  values
}

aaa_marker_hits <- function(values, genes) {
  sum(values[genes] %in% TRUE)
}

aaa_marker_all <- function(values, genes) {
  length(genes) > 0L && all(values[genes] %in% TRUE)
}

aaa_marker_any <- function(values, genes) {
  any(values[genes] %in% TRUE)
}


#' Classify canonical proline biosynthesis
classify_proline_biosynthesis <- function(x) {
  genes <- c("proA", "proB", "proC")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown proline biosynthetic potential")
  }
  hits <- aaa_marker_hits(values, genes)
  if (hits == 3L) {
    return("Complete proline biosynthetic potential")
  }
  if (hits == 2L) {
    return("Probable proline biosynthetic potential")
  }
  if (hits == 1L) {
    return("Partial proline biosynthetic evidence")
  }
  "No detected proline biosynthetic potential"
}


#' Classify pipecolate-associated biosynthesis
#'
#' Pipecolate can be produced by several enzyme routes. A strong call therefore
#' requires coherent evidence from at least two route steps and is deliberately
#' described as potential rather than as a confirmed complete pathway.
classify_pipecolate_biosynthesis <- function(x) {
  genes <- c("lhpI", "dpkA", "lysDH", "proC", "p2cr")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown pipecolate biosynthetic potential")
  }

  route_entry <- aaa_marker_any(values, c("lhpI", "dpkA", "lysDH"))
  route_reduction <- aaa_marker_any(values, c("proC", "p2cr"))
  hits <- aaa_marker_hits(values, genes)

  if (route_entry && route_reduction && hits >= 2L) {
    return("Pipecolate biosynthetic potential")
  }
  if (hits >= 2L) {
    return("Probable pipecolate-associated pathway potential")
  }
  if (hits == 1L) {
    return("Partial pipecolate-pathway evidence")
  }
  "No detected pipecolate biosynthetic potential"
}


#' Classify carotenoid pigment-production potential
#'
#' The classifier is hierarchical. A downstream pigment is only reported when
#' its precursor route and product-specific conversion markers are present.
classify_carotenoid_pigments <- function(x) {
  genes <- c(
    "crtE", "crtB", "crtP", "zIso", "crtQ", "crtH", "crtI", "al1",
    "crtY", "cruA", "cruP", "al2", "crtR", "crtZ", "lut5",
    "zep", "nsy", "crtW", "crtC", "crtD", "crtF",
    "crtLe", "crtLb", "lut1", "crtO", "carT", "carD",
    "lyeJ", "cruF"
  )
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown carotenoid pigment potential")
  }

  precursor <- aaa_marker_all(values, c("crtE", "crtB"))
  lycopene_conversion <- (
    aaa_marker_any(values, c("crtI", "al1")) ||
      aaa_marker_all(values, c("crtP", "zIso", "crtQ", "crtH"))
  )
  lycopene <- precursor && lycopene_conversion

  cyclase <- aaa_marker_any(values, c("crtY", "cruA", "cruP", "al2"))
  beta_carotene <- lycopene && cyclase
  hydroxylase <- aaa_marker_any(values, c("crtZ", "crtR", "lut5"))
  zeaxanthin <- beta_carotene && hydroxylase

  products <- character()
  if (lycopene) products <- c(products, "lycopene")
  if (beta_carotene) products <- c(products, "beta-carotene")
  if (zeaxanthin) products <- c(products, "zeaxanthin")

  if (zeaxanthin && isTRUE(values[["zep"]])) {
    products <- c(products, "violaxanthin")
  }
  if (zeaxanthin && isTRUE(values[["zep"]]) && isTRUE(values[["nsy"]])) {
    products <- c(products, "neoxanthin")
  }
  if (beta_carotene && hydroxylase && isTRUE(values[["crtW"]])) {
    products <- c(products, "astaxanthin")
  }
  if (beta_carotene && aaa_marker_any(values, c("crtW", "crtO"))) {
    products <- c(products, "canthaxanthin")
  }
  if (lycopene && aaa_marker_all(values, c("crtC", "crtD", "crtF"))) {
    products <- c(products, "spirilloxanthin")
  }
  if (lycopene &&
    aaa_marker_all(values, c("crtLe", "crtLb", "lut5", "lut1"))) {
    products <- c(products, "lutein")
  }
  if (lycopene && aaa_marker_all(values, c("al1", "al2", "carT", "carD"))) {
    products <- c(products, "neurosporaxanthin")
  }
  if (lycopene && aaa_marker_all(values, c("lyeJ", "crtD", "cruF"))) {
    products <- c(products, "bacterioruberin")
  }

  if (length(products)) {
    return(paste(
      "Carotenoid pigment potential:",
      paste(unique(products), collapse = "; ")
    ))
  }

  hits <- aaa_marker_hits(values, genes)
  if (hits >= 3L) {
    return("Probable carotenoid-pathway potential")
  }
  if (hits >= 1L) {
    return("Partial carotenoid-pathway evidence")
  }
  "No detected carotenoid pigment potential"
}


#' Classify bacteriorhodopsin-based phototrophy
classify_bacteriorhodopsin <- function(x) {
  genes <- c("bop", "blh", "bcmo1")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown bacteriorhodopsin potential")
  }

  opsin <- isTRUE(values[["bop"]])
  retinal <- aaa_marker_any(values, c("blh", "bcmo1"))

  if (opsin && retinal) {
    return("Bacteriorhodopsin phototrophy potential")
  }
  if (opsin) {
    return("Bacteriorhodopsin protein without confirmed retinal-production marker")
  }
  if (retinal) {
    return("Retinal-production potential without bacteriorhodopsin marker")
  }
  "No detected bacteriorhodopsin potential"
}


#' Classify syngas-to-ethanol valorisation potential
#'
#' This is distinct from the existing homoacetogenesis and CO-to-H2 entries:
#' it requires a coherent Wood-Ljungdahl/CO-utilisation backbone together with
#' an ethanol-forming branch.
classify_syngas_to_ethanol <- function(x) {
  genes <- c(
    "cooS", "cooF", "acsA", "acsB", "acsC", "acsD", "acsE",
    "fhs", "folD", "metF", "fdh", "adhE", "aor", "adh"
  )
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown syngas-to-ethanol potential")
  }

  co_entry <- aaa_marker_any(values, c("cooS", "acsA"))
  carbonyl <- isTRUE(values[["acsB"]]) &&
    aaa_marker_hits(values, c("acsC", "acsD", "acsE")) >= 2L
  methyl <- aaa_marker_hits(values, c("fhs", "folD", "metF", "fdh")) >= 3L
  ethanol_branch <- isTRUE(values[["adhE"]]) ||
    aaa_marker_all(values, c("aor", "adh"))

  if (co_entry && carbonyl && methyl && ethanol_branch) {
    return("Syngas-to-ethanol valorisation potential")
  }
  if (carbonyl && methyl && ethanol_branch) {
    return("Probable C1-to-ethanol valorisation potential")
  }
  if (ethanol_branch &&
    (co_entry || carbonyl || methyl)) {
    return("Partial syngas-to-ethanol pathway evidence")
  }
  if (co_entry || carbonyl || methyl) {
    return("C1 assimilation evidence without a resolved ethanol branch")
  }
  "No detected syngas-to-ethanol potential"
}


#' Classify reductive glycine pathway potential
#'
#' This pathway supports assimilation of formate/CO2 into glycine, serine and
#' ultimately pyruvate. The call requires the formate-THF module, reverse
#' glycine-cleavage machinery, serine formation and a serine-to-pyruvate step.
classify_reductive_glycine_pathway <- function(x) {
  genes <- c(
    "fdh", "fhs", "folD",
    "gcvP", "gcvT", "gcvH", "lpd",
    "glyA", "sdaA", "sdaB", "tdcG"
  )
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown reductive glycine pathway potential")
  }

  formate_module <- aaa_marker_hits(values, c("fdh", "fhs", "folD")) >= 2L
  reverse_gcs <- aaa_marker_hits(values, c("gcvP", "gcvT", "gcvH", "lpd")) >= 3L
  serine_module <- isTRUE(values[["glyA"]])
  pyruvate_module <- aaa_marker_any(values, c("sdaA", "sdaB", "tdcG"))

  if (formate_module && reverse_gcs && serine_module && pyruvate_module) {
    return("Complete reductive glycine pathway potential")
  }
  if (formate_module && reverse_gcs && serine_module) {
    return("Reductive glycine pathway potential to serine")
  }

  modules <- sum(c(formate_module, reverse_gcs, serine_module, pyruvate_module))
  if (modules >= 2L) {
    return("Probable reductive glycine pathway potential")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Partial reductive glycine pathway evidence")
  }
  "No detected reductive glycine pathway potential"
}

# =============================================================================
# Additional methane/C1 utilisation, dark fermentation and osmoprotection
# classifiers. These entries complement the existing methanotrophy, ectoine,
# proline, pipecolate, CO and syngas functions without replacing them.
# =============================================================================

#' Classify ribulose-monophosphate formaldehyde assimilation
#'
#' Hps and Phi catalyse the two defining formaldehyde-fixation reactions of
#' the RuMP pathway. Methane oxidation is intentionally not inferred here:
#' that is handled by the pre-existing methanotrophy classifier.
classify_rump_formaldehyde_assimilation <- function(x) {
  genes <- c("hps", "phi", "fae")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown RuMP formaldehyde-assimilation potential")
  }

  core <- aaa_marker_all(values, c("hps", "phi"))
  if (core && isTRUE(values[["fae"]])) {
    return("Complete RuMP formaldehyde-assimilation potential with formaldehyde-handling support")
  }
  if (core) {
    return("Complete RuMP formaldehyde-assimilation potential")
  }
  if (aaa_marker_hits(values, c("hps", "phi")) == 1L) {
    return("Partial RuMP formaldehyde-assimilation evidence")
  }
  if (isTRUE(values[["fae"]])) {
    return("Formaldehyde-handling evidence without a complete RuMP core")
  }
  "No detected RuMP formaldehyde-assimilation potential"
}


#' Classify formate-hydrogen-lyase dark-fermentation potential
#'
#' A strong call requires the formate dehydrogenase-H component and the
#' hydrogenase-3 catalytic/electron-transfer core. Regulatory and membrane
#' subunits strengthen, but do not independently establish, the pathway.
classify_dark_fermentation_fhl <- function(x) {
  genes <- c(
    "fdhF", "hycE", "hycB", "hycF", "hycG",
    "hycC", "hycD", "hycI", "fhlA"
  )
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown formate-hydrogen-lyase potential")
  }

  donor <- isTRUE(values[["fdhF"]])
  catalytic <- isTRUE(values[["hycE"]])
  transfer <- aaa_marker_hits(values, c("hycB", "hycF", "hycG")) >= 2L
  assembly <- aaa_marker_any(values, c("hycI", "fhlA"))
  membrane <- aaa_marker_hits(values, c("hycC", "hycD")) >= 1L

  if (donor && catalytic && transfer && assembly) {
    return("High-confidence formate-hydrogen-lyase dark-fermentation potential")
  }
  if (donor && catalytic && transfer) {
    return("Formate-hydrogen-lyase dark-fermentation potential")
  }

  modules <- sum(c(donor, catalytic, transfer, assembly, membrane))
  if (modules >= 3L) {
    return("Probable formate-hydrogen-lyase potential")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Partial formate-hydrogen-lyase evidence")
  }
  "No detected formate-hydrogen-lyase potential"
}


#' Classify PFOR/ferredoxin-hydrogenase dark-fermentation potential
#'
#' The classifier requires a pyruvate:ferredoxin oxidoreductase route and a
#' ferredoxin-linked [FeFe]-hydrogenase module. Either a monomeric NifJ-type
#' PFOR or a sufficiently complete PorABDG complex can satisfy the PFOR arm.
classify_dark_fermentation_pfor_hydrogenase <- function(x) {
  genes <- c(
    "nifJ", "porA", "porB", "porD", "porG",
    "hydA", "hydB", "hydC", "fdx"
  )
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown PFOR-linked dark-fermentation potential")
  }

  pfor <- isTRUE(values[["nifJ"]]) ||
    aaa_marker_hits(values, c("porA", "porB", "porD", "porG")) >= 3L
  hydrogenase <- isTRUE(values[["hydA"]]) &&
    aaa_marker_hits(values, c("hydB", "hydC")) >= 1L
  electron_carrier <- isTRUE(values[["fdx"]])

  if (pfor && hydrogenase && electron_carrier) {
    return("High-confidence PFOR-linked dark-fermentative hydrogen potential")
  }
  if (pfor && hydrogenase) {
    return("PFOR-linked dark-fermentative hydrogen potential")
  }
  if (pfor || hydrogenase) {
    return("Partial PFOR/hydrogenase dark-fermentation evidence")
  }
  "No detected PFOR-linked dark-fermentation potential"
}


#' Classify glycine-betaine biosynthesis from choline
#'
#' The route is called complete when both choline-to-betaine-aldehyde and
#' betaine-aldehyde-to-glycine-betaine activities are represented. Alternative
#' bacterial naming conventions (BetA/BetB and GbsB/GbsA) are accepted.
classify_glycine_betaine_biosynthesis <- function(x) {
  genes <- c("betA", "betB", "gbsA", "gbsB")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown glycine-betaine biosynthetic potential")
  }

  first_step <- aaa_marker_any(values, c("betA", "gbsB"))
  second_step <- aaa_marker_any(values, c("betB", "gbsA"))

  if (first_step && second_step) {
    return("Complete glycine-betaine biosynthetic potential from choline")
  }
  if (first_step || second_step) {
    return("Partial glycine-betaine biosynthetic evidence")
  }
  "No detected glycine-betaine biosynthetic potential"
}


#' Classify trehalose biosynthesis
#'
#' OtsAB and TreYZ are treated as independent complete routes. TreS alone is
#' reported as an alternative trehalose-interconversion potential rather than
#' as de novo biosynthesis.
classify_trehalose_biosynthesis <- function(x) {
  genes <- c("otsA", "otsB", "treY", "treZ", "treS")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown trehalose biosynthetic potential")
  }

  otsab <- aaa_marker_all(values, c("otsA", "otsB"))
  treyz <- aaa_marker_all(values, c("treY", "treZ"))
  tres <- isTRUE(values[["treS"]])

  if (otsab && treyz) {
    return("Multiple complete trehalose biosynthetic routes")
  }
  if (otsab) {
    return("Complete OtsAB trehalose biosynthetic potential")
  }
  if (treyz) {
    return("Complete TreYZ trehalose biosynthetic potential")
  }
  if (tres) {
    return("TreS-dependent trehalose interconversion potential")
  }

  if (aaa_marker_hits(values, c("otsA", "otsB", "treY", "treZ")) > 0L) {
    return("Partial trehalose biosynthetic evidence")
  }
  "No detected trehalose biosynthetic potential"
}


#' Classify mannosylglycerate biosynthesis
#'
#' The two-step MpgS/MpgP route and the bifunctional MpgS/P route are accepted.
#' The metabolite is interpreted as a compatible-solute potential, not as proof
#' of accumulation under osmotic or thermal stress.
classify_mannosylglycerate_biosynthesis <- function(x) {
  genes <- c("mpgS", "mpgP", "mpgSP")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown mannosylglycerate biosynthetic potential")
  }

  if (isTRUE(values[["mpgSP"]])) {
    return("Complete bifunctional mannosylglycerate biosynthetic potential")
  }
  if (aaa_marker_all(values, c("mpgS", "mpgP"))) {
    return("Complete two-step mannosylglycerate biosynthetic potential")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Partial mannosylglycerate biosynthetic evidence")
  }
  "No detected mannosylglycerate biosynthetic potential"
}

# =============================================================================
# Lactate, lactate-derived biopolymers, C1 assimilation and additional
# compatible-solute classifiers.
# =============================================================================

#' Classify L-lactate production
#'
#' A stereospecific L-lactate dehydrogenase is required. Generic ldh
#' annotations are intentionally not used as diagnostic evidence because they
#' may denote respiratory or non-stereospecific enzymes.
classify_l_lactate_production <- function(x) {
  genes <- c("ldhL")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown L-lactate production potential")
  }
  if (isTRUE(values[["ldhL"]])) {
    return("L-lactate production potential")
  }
  "No detected L-lactate production potential"
}


#' Classify D-lactate production
classify_d_lactate_production <- function(x) {
  genes <- c("ldhD")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown D-lactate production potential")
  }
  if (isTRUE(values[["ldhD"]])) {
    return("D-lactate production potential")
  }
  "No detected D-lactate production potential"
}


#' Classify respiratory lactate utilisation
#'
#' Complete LutABC/LldEFG complexes are treated as strong evidence. Standalone
#' stereospecific respiratory dehydrogenases are reported as narrower L- or
#' D-lactate oxidation potential and are not promoted to a complete complex.
classify_respiratory_lactate_utilisation <- function(x) {
  genes <- c(
    "lutA", "lutB", "lutC", "lutP",
    "lldE", "lldF", "lldG", "lldP", "lldD", "dld"
  )
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown respiratory lactate-utilisation potential")
  }

  lut <- aaa_marker_all(values, c("lutA", "lutB", "lutC"))
  lldefg <- aaa_marker_all(values, c("lldE", "lldF", "lldG"))
  l_specific <- isTRUE(values[["lldD"]])
  d_specific <- isTRUE(values[["dld"]])
  transporter <- aaa_marker_any(values, c("lutP", "lldP"))

  if ((lut || lldefg) && transporter) {
    return("High-confidence respiratory lactate-utilisation potential")
  }
  if (lut || lldefg) {
    return("Respiratory lactate-utilisation potential")
  }
  if (l_specific && d_specific) {
    return("L- and D-lactate oxidation potential")
  }
  if (l_specific) {
    return("L-lactate oxidation potential")
  }
  if (d_specific) {
    return("D-lactate oxidation potential")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Partial lactate-utilisation evidence")
  }
  "No detected respiratory lactate-utilisation potential"
}


#' Classify conversion of lactate carbon into PHA
#'
#' This predicts a metabolic chassis able to oxidise lactate to central carbon
#' and carry a canonical PHB/PHA module. It does not imply direct incorporation
#' of lactyl-CoA into the polymer.
classify_lactate_to_pha <- function(x) {
  genes <- c(
    "lutA", "lutB", "lutC", "lldE", "lldF", "lldG", "lldD", "dld",
    "phaA", "phaB", "phaC"
  )
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown lactate-to-PHA potential")
  }

  lactate <- aaa_marker_all(values, c("lutA", "lutB", "lutC")) ||
    aaa_marker_all(values, c("lldE", "lldF", "lldG")) ||
    aaa_marker_any(values, c("lldD", "dld"))
  pha_core <- aaa_marker_all(values, c("phaA", "phaB", "phaC"))

  if (lactate && pha_core) {
    return("Lactate-fed PHA biosynthetic potential")
  }
  if (lactate && aaa_marker_hits(values, c("phaA", "phaB", "phaC")) >= 2L) {
    return("Probable lactate-fed PHA biosynthetic potential")
  }
  if (lactate || aaa_marker_hits(values, c("phaA", "phaB", "phaC")) > 0L) {
    return("Partial lactate-to-PHA pathway evidence")
  }
  "No detected lactate-to-PHA potential"
}


#' Classify a putative PLA/lactate-copolyester engineering chassis
#'
#' Natural phaC annotations do not establish lactyl-CoA polymerisation. The
#' classifier therefore reports an engineering chassis only and requires a
#' lactate-forming enzyme, propionyl-CoA transferase and PHA synthase.
classify_pla_engineering_chassis <- function(x) {
  genes <- c("ldhL", "ldhD", "pct", "phaC")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown PLA engineering-chassis potential")
  }

  lactate_supply <- aaa_marker_any(values, c("ldhL", "ldhD"))
  lactyl_coa_step <- isTRUE(values[["pct"]])
  polymerase <- isTRUE(values[["phaC"]])

  if (lactate_supply && lactyl_coa_step && polymerase) {
    return("Putative PLA/lactate-copolyester engineering chassis")
  }
  if (sum(c(lactate_supply, lactyl_coa_step, polymerase)) == 2L) {
    return("Partial PLA engineering-chassis evidence")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Limited PLA engineering-chassis evidence")
  }
  "No detected PLA engineering-chassis potential"
}


#' Classify methanol oxidation
#'
#' PQQ-dependent MxaFI and lanthanide-dependent XoxF systems are evaluated
#' separately. Accessory cytochrome and assembly proteins strengthen the call.
classify_methanol_oxidation <- function(x) {
  genes <- c("mxaF", "mxaI", "mxaG", "mxaJ", "xoxF", "xoxG", "xoxJ")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown methanol-oxidation potential")
  }

  mxa_core <- aaa_marker_all(values, c("mxaF", "mxaI"))
  mxa_support <- aaa_marker_any(values, c("mxaG", "mxaJ"))
  xox_core <- isTRUE(values[["xoxF"]])
  xox_support <- aaa_marker_any(values, c("xoxG", "xoxJ"))

  if (mxa_core && mxa_support && xox_core && xox_support) {
    return("MxaFI and XoxF methanol-oxidation potential")
  }
  if (mxa_core && mxa_support) {
    return("High-confidence MxaFI methanol-oxidation potential")
  }
  if (xox_core && xox_support) {
    return("High-confidence XoxF methanol-oxidation potential")
  }
  if (mxa_core) {
    return("MxaFI methanol-oxidation potential")
  }
  if (xox_core) {
    return("XoxF-associated methanol-oxidation potential")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Partial methanol-oxidation evidence")
  }
  "No detected methanol-oxidation potential"
}


#' Classify serine-cycle C1 assimilation
#'
#' A strong call requires formaldehyde incorporation, glyoxylate regeneration
#' and the characteristic malyl-CoA cleavage arm. Common central enzymes such
#' as GlyA are supporting rather than independently diagnostic.
classify_serine_cycle_c1_assimilation <- function(x) {
  genes <- c("fae", "mtdA", "mch", "glyA", "hprA", "sgaA", "mtkA", "mtkB", "mcl")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown serine-cycle C1-assimilation potential")
  }

  c1_entry <- aaa_marker_hits(values, c("fae", "mtdA", "mch", "glyA")) >= 3L
  serine_to_glycerate <- aaa_marker_all(values, c("hprA", "sgaA"))
  malyl_coa_arm <- isTRUE(values[["mcl"]]) && aaa_marker_hits(values, c("mtkA", "mtkB")) >= 1L

  if (c1_entry && serine_to_glycerate && malyl_coa_arm) {
    return("Complete serine-cycle C1-assimilation potential")
  }
  if (sum(c(c1_entry, serine_to_glycerate, malyl_coa_arm)) >= 2L) {
    return("Probable serine-cycle C1-assimilation potential")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Partial serine-cycle evidence")
  }
  "No detected serine-cycle C1-assimilation potential"
}


#' Classify glucosylglycerol biosynthesis
classify_glucosylglycerol_biosynthesis <- function(x) {
  genes <- c("ggpS", "ggpP", "ggpSP")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown glucosylglycerol biosynthetic potential")
  }
  if (isTRUE(values[["ggpSP"]])) {
    return("Complete bifunctional glucosylglycerol biosynthetic potential")
  }
  if (aaa_marker_all(values, c("ggpS", "ggpP"))) {
    return("Complete glucosylglycerol biosynthetic potential")
  }
  if (aaa_marker_hits(values, genes) > 0L) {
    return("Partial glucosylglycerol biosynthetic evidence")
  }
  "No detected glucosylglycerol biosynthetic potential"
}


#' Classify N-acetylglutaminylglutamine amide biosynthesis
classify_naggn_biosynthesis <- function(x) {
  genes <- c("asnO", "ngg")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown NAGGN biosynthetic potential")
  }
  if (aaa_marker_all(values, genes)) {
    return("Complete NAGGN biosynthetic potential")
  }
  if (aaa_marker_hits(values, genes) == 1L) {
    return("Partial NAGGN biosynthetic evidence")
  }
  "No detected NAGGN biosynthetic potential"
}


#' Classify sucrose biosynthesis as a compatible-solute route
classify_sucrose_biosynthesis <- function(x) {
  genes <- c("sps", "spp")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown sucrose biosynthetic potential")
  }
  if (aaa_marker_all(values, genes)) {
    return("Complete SPS/SPP sucrose biosynthetic potential")
  }
  if (aaa_marker_hits(values, genes) == 1L) {
    return("Partial sucrose biosynthetic evidence")
  }
  "No detected sucrose biosynthetic potential"
}

# =============================================================================
# Iron-cycle classifiers.
#
# Dissimilatory iron reduction and iron oxidation are deliberately conservative:
# the diagnostic cytochromes are distributed across different lineages
# (Shewanella omcA, Geobacter omcS/omcZ; Fe(II) oxidiser cyc2), so no single
# genome is expected to carry the whole set. The classifiers therefore never
# call a "Complete" pathway from these markers, only "Probable"/"Partial",
# which avoids over-interpreting individual outer-membrane cytochrome hits.
# Genes prone to name collisions across pathways (the methanogen mtrA-H
# methyltransferase shares symbols with metal-reduction Mtr) are excluded.
# =============================================================================

#' Classify dissimilatory (respiratory) iron reduction
classify_dissimilatory_iron_reduction <- function(x) {
  genes <- c("omcA", "omcB", "omcS", "omcZ")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown dissimilatory iron-reduction potential")
  }
  hits <- aaa_marker_hits(values, genes)
  if (hits >= 2L) {
    return("Probable dissimilatory iron-reduction potential")
  }
  if (hits == 1L) {
    return("Partial dissimilatory iron-reduction evidence")
  }
  "No detected dissimilatory iron-reduction potential"
}

#' Classify iron (Fe(II)) oxidation
classify_iron_oxidation <- function(x) {
  genes <- c("cyc2", "cyc1")
  values <- aaa_marker_vector(x, genes)
  if (all(is.na(values))) {
    return("Unknown iron-oxidation potential")
  }
  if (values["cyc2"] %in% TRUE) {
    return("Probable iron-oxidation potential")
  }
  if (aaa_marker_any(values, genes)) {
    return("Partial iron-oxidation evidence")
  }
  "No detected iron-oxidation potential"
}
