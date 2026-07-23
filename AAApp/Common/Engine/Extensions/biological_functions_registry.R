# =============================================================================
# Central registry of biological functions
#
# Each function defines:
#   - genes searched in GFF annotations;
#   - curated roles: diagnostic, supporting and accessory;
#   - classifier;
#   - functional-abundance selector;
#   - descriptions shown automatically in Shiny and documentation.
#
# Uses NCBI Taxonomy, RefSeq assemblies and GFF annotations.
# Biological interpretation is determined exclusively by the curated registry.
# =============================================================================

aaa_gene_roles <- function(diagnostic = character(),
                           supporting = character(),
                           accessory = character()) {
  list(
    diagnostic = unique(diagnostic),
    supporting = unique(supporting),
    accessory = unique(accessory)
  )
}

aaa_registry_genes <- function(gene_roles) {
  unique(unlist(gene_roles, use.names = FALSE))
}

biological_function_registry <- list(
  Methanogenesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("mcrA", "mcrB", "mcrG")
    )
    list(
      display_name = "Methanogenesis (general)",
      category = "Methane and C1 metabolism",
      description = paste(
        "Detects the shared terminal methyl-coenzyme M reductase complex.",
        "It identifies general methanogenic potential but does not resolve",
        "the substrate route."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_methanogenesis,
      graph_main = "General methanogenic potential",
      analysis_name = "Methanogenesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "mcr genes are diagnostic of methane formation but are shared by",
        "major methanogenic routes. They do not distinguish hydrogenotrophic",
        "from acetoclastic methanogenesis."
      )
    )
  }),
  Methanogenesis_subtype = local({
    roles <- aaa_gene_roles(
      diagnostic = c(
        "mcrA", "mcrB", "mcrG",
        "ftr", "mch", "mtd", "mer",
        "ackA", "pta", "acs",
        "cdhA", "cdhB", "cdhC", "cdhD", "cdhE",
        "mtaA", "mtaB", "mtaC",
        "mtbA", "mtmB", "mtmC", "mtbB", "mtbC", "mttB", "mttC",
        "acsA", "acsB", "acsC", "acsD", "acsE",
        "fhs", "folD", "metF", "fdh"
      ),
      supporting = c(
        "fwdA", "fwdB", "fwdC", "fwdD",
        "fmdA", "fmdB", "fmdC", "fmdD",
        "mtrA", "mtrB", "mtrC", "mtrD",
        "mtrE", "mtrF", "mtrG", "mtrH"
      ),
      accessory = character()
    )
    list(
      display_name = "Methanogenesis and Wood-Ljungdahl profile",
      category = "Methane and C1 metabolism",
      description = paste(
        "Classifies complete hydrogenotrophic, acetoclastic and methylotrophic",
        "methanogenesis routes. A complete reductive Wood-Ljungdahl profile is",
        "reported separately only when methanogenesis is not supported; partial",
        "profiles are grouped as insufficient evidence."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      diagnostic_groups = list(
        Shared_methanogenic_terminal_step = c("mcrA", "mcrB", "mcrG"),
        Hydrogenotrophic_CO2_reduction = c(
          "fwdA", "fwdB", "fwdC", "fwdD",
          "fmdA", "fmdB", "fmdC", "fmdD",
          "ftr", "mch", "mtd", "mer"
        ),
        Methyl_transfer = c(
          "mtrA", "mtrB", "mtrC", "mtrD",
          "mtrE", "mtrF", "mtrG", "mtrH"
        ),
        Acetoclastic_route = c(
          "ackA", "pta", "acs",
          "cdhA", "cdhB", "cdhC", "cdhD", "cdhE"
        ),
        Methylotrophic_substrate_entry = c(
          "mtaA", "mtaB", "mtaC",
          "mtbA", "mtmB", "mtmC", "mtbB", "mtbC", "mttB", "mttC"
        ),
        Homoacetogenic_Wood_Ljungdahl = c(
          "acsA", "acsB", "acsC", "acsD", "acsE",
          "fhs", "folD", "metF", "fdh"
        )
      ),
      classifier = classify_methanogenesis_subtype,
      graph_main = "Methanogenesis and Wood-Ljungdahl profile",
      graph_note = paste(
        "A complete reductive Wood-Ljungdahl pathway detected without",
        "methanogenesis-specific evidence may indicate homoacetogenic potential."
      ),
      analysis_name = "Methanogenesis_subtype",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Methanogenesis takes precedence over homoacetogenesis when both Mcr",
        "markers and Wood-Ljungdahl genes are detected. Methylotrophic",
        "methanogenesis requires a supported Mcr complex and a coherent methanol",
        "or methylamine entry module. Partial routes are classified as",
        "insufficient evidence."
      )
    )
  }),
  Hydrogenotrophic_methanogenesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("mcrA", "mcrB", "mcrG", "ftr", "mch", "mtd", "mer"),
      supporting = c(
        "fwdA", "fwdB", "fwdC", "fwdD",
        "fmdA", "fmdB", "fmdC", "fmdD",
        "mtrA", "mtrB", "mtrC", "mtrD",
        "mtrE", "mtrF", "mtrG", "mtrH"
      )
    )
    list(
      display_name = "Hydrogenotrophic methanogenesis",
      category = "Methane and C1 metabolism",
      description = paste(
        "Evaluates genomic potential for CO2-reducing hydrogenotrophic",
        "methanogenesis using the shared mcr terminal complex together",
        "with CO2-reduction and methyl-transfer pathway evidence."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_hydrogenotrophic_methanogenesis,
      graph_main = "Hydrogenotrophic methanogenesis potential",
      analysis_name = "Hydrogenotrophic_methanogenesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Genomic potential does not demonstrate in-situ activity."
    )
  }),
  Acetoclastic_methanogenesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c(
        "mcrA", "mcrB", "mcrG", "ackA", "pta", "acs",
        "cdhA", "cdhB", "cdhC", "cdhD", "cdhE"
      )
    )
    list(
      display_name = "Acetoclastic methanogenesis",
      category = "Methane and C1 metabolism",
      description = paste(
        "Evaluates genomic potential for acetate-dependent methanogenesis",
        "using the shared mcr terminal complex, acetate activation and",
        "CODH/ACS complex evidence."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_acetoclastic_methanogenesis,
      graph_main = "Acetoclastic methanogenesis potential",
      analysis_name = "Acetoclastic_methanogenesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Acetate-activation genes alone are not sufficient; the classifier",
        "also requires methanogenesis-specific pathway evidence."
      )
    )
  }),
  Homoacetogenesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("acsA", "acsB", "acsC", "acsD"),
      supporting = c("acsE", "fhs", "folD", "metF", "fdh"),
      accessory = c("pta", "ackA", "mcrA", "mcrB", "mcrG")
    )
    list(
      display_name = "Homoacetogenesis (Wood-Ljungdahl pathway)",
      category = "Carbon fixation and CO/CO2",
      description = paste(
        "Evaluates a complete reductive Wood-Ljungdahl profile as possible",
        "homoacetogenesis only when a supported methanogenic Mcr complex is",
        "absent. Partial profiles are reported as insufficient evidence."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      diagnostic_groups = list(
        Carbonyl_branch = c("acsA", "acsB"),
        Corrinoid_iron_sulfur_module = c("acsC", "acsD", "acsE"),
        Methyl_branch = c("fhs", "folD", "metF"),
        Electron_entry = c("fdh"),
        Acetate_output = c("pta", "ackA"),
        Methanogenesis_exclusion = c("mcrA", "mcrB", "mcrG")
      ),
      classifier = classify_homoacetogenesis,
      evidence_scorer = score_homoacetogenesis,
      graph_main = "Homoacetogenic potential",
      analysis_name = "Homoacetogenesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "A complete Wood-Ljungdahl pathway is not exclusive to homoacetogens.",
        "When a supported Mcr terminal complex is present, homoacetogenesis is",
        "not assigned. The result represents genomic potential, not activity."
      )
    )
  }),
  Methanotrophy = local({
    roles <- aaa_gene_roles(
      diagnostic = c(
        "pmoA", "pmoB", "pmoC",
        "mmoX", "mmoY", "mmoZ"
      ),
      supporting = c("mmoB", "mmoC", "mmoD")
    )
    list(
      display_name = "Methanotrophy",
      category = "Methane and C1 metabolism",
      description = paste(
        "Detects particulate and soluble methane monooxygenase systems.",
        "The analysis distinguishes complete, probable and partial evidence."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_methanotroph,
      graph_main = "Methanotrophic potential",
      analysis_name = "Methanotrophy",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Methane monooxygenase genes provide the principal diagnostic",
        "evidence used by this classifier."
      )
    )
  }),
  Ectoine = local({
    roles <- aaa_gene_roles(
      diagnostic = c("ectA", "ectB", "ectC"),
      supporting = "ectD"
    )
    list(
      display_name = "Ectoine production",
      category = "Osmoprotectors and high value molecules",
      description = paste(
        "Detects the ectABC biosynthetic pathway and ectD-dependent",
        "hydroxyectoine formation."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_ectoine,
      graph_main = "Ectoine-production potential",
      analysis_name = "Ectoine",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Nitrification = local({
    roles <- aaa_gene_roles(
      diagnostic = c("amoA", "amoB", "amoC", "nxrA", "nxrB"),
      supporting = "hao"
    )
    list(
      display_name = "Nitrification",
      category = "Nitrogen cycle",
      description = paste(
        "Evaluates ammonia oxidation and nitrite oxidation markers.",
        "amo genes are analyzed as nitrification markers, not",
        "methanotrophy markers."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_nitrification,
      graph_main = "Nitrification potential",
      analysis_name = "Nitrification",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Denitrification = local({
    roles <- aaa_gene_roles(
      diagnostic = c("narG", "nirK", "nirS", "norB", "nosZ")
    )
    list(
      display_name = "Denitrification",
      category = "Nitrogen cycle",
      description = paste(
        "Evaluates nitrate, nitrite, nitric oxide and nitrous oxide",
        "reduction markers."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_denitrification,
      graph_main = "Denitrification potential",
      analysis_name = "Denitrification",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Sulfate_reduction = local({
    roles <- aaa_gene_roles(
      diagnostic = c("aprA", "aprB", "dsrA", "dsrB"),
      supporting = "sat"
    )
    list(
      display_name = "Dissimilatory sulfate reduction",
      category = "Sulfur cycle",
      description = paste(
        "Evaluates sulfate activation, APS reduction and dissimilatory",
        "sulfite reduction markers."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_sulfate_reducer,
      graph_main = "Sulfate-reduction potential",
      analysis_name = "Sulfate_reduction",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Nitrogen_fixation = local({
    roles <- aaa_gene_roles(
      diagnostic = c("nifH", "nifD", "nifK")
    )
    list(
      display_name = "Nitrogen fixation",
      category = "Nitrogen cycle",
      description = "Evaluates the nifHDK nitrogenase core.",
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_nitrogen_fixation,
      graph_main = "Nitrogen-fixation potential",
      analysis_name = "Nitrogen_fixation",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Anammox = local({
    roles <- aaa_gene_roles(
      diagnostic = c("hzsA", "hzsB", "hzsC", "hdh")
    )
    list(
      display_name = "Anammox",
      category = "Nitrogen cycle",
      description = paste(
        "Evaluates hydrazine synthase and hydrazine dehydrogenase",
        "markers associated with anaerobic ammonium oxidation."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_anammox,
      graph_main = "Anammox potential",
      analysis_name = "Anammox",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Acetogenesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("acsB", "fhs"),
      supporting = c("acsA", "fdh")
    )
    list(
      display_name = "Acetogenesis",
      category = "Carbon fixation and CO/CO2",
      description = paste(
        "Evaluates selected Wood-Ljungdahl pathway markers.",
        "This entry is distinct from acetoclastic methanogenesis."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_acetogenesis,
      graph_main = "Acetogenic potential",
      analysis_name = "Acetogenesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "The marker set is intentionally restricted to improve",
        "the specificity of the inferred acetogenic potential."
      )
    )
  }),
  DNRA = local({
    roles <- aaa_gene_roles(diagnostic = c("nrfA", "nrfH"), supporting = c("narG", "napA"))
    list(
      display_name = "DNRA (nitrate/nitrite to ammonium)", category = "Nitrogen cycle",
      description = "Evaluates dissimilatory nitrate reduction to ammonium using the nrfAH nitrite-reductase module.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("DNRA", list(Nitrite_to_ammonium = c("nrfA", "nrfH")), list(Nitrate_entry = c("narG", "napA")), .5, .5),
      graph_main = "DNRA potential", analysis_name = "DNRA",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "nrfA is the main diagnostic marker; nitrate-reductase genes alone are not sufficient."
    )
  }),
  Nitrous_oxide_reduction = local({
    roles <- aaa_gene_roles(diagnostic = c("nosZ"), supporting = c("nosD", "nosF", "nosY", "nosL"))
    list(
      display_name = "Nitrous oxide reduction", category = "Nitrogen cycle",
      description = "Evaluates genomic potential for the final N2O-to-N2 step of denitrification.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("nitrous oxide reduction", list(Catalytic = c("nosZ")), list(Maturation = c("nosD", "nosF", "nosY", "nosL")), 1, 1),
      graph_main = "N2O-reduction potential", analysis_name = "Nitrous_oxide_reduction",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Polyphosphate_metabolism = local({
    roles <- aaa_gene_roles(diagnostic = c("ppk1", "ppk2"), supporting = c("ppx", "phoU"))
    list(
      display_name = "Polyphosphate accumulation and turnover", category = "Phosphorus cycle",
      description = "Evaluates polyphosphate kinase and exopolyphosphatase markers relevant to EBPR, without assigning a PAO phenotype from taxonomy alone.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("polyphosphate metabolism", list(PolyP_synthesis = c("ppk1", "ppk2")), list(PolyP_turnover = c("ppx"), Regulation = c("phoU")), .5, .5),
      graph_main = "Polyphosphate-metabolism potential", analysis_name = "Polyphosphate_metabolism",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "These markers support polyphosphate metabolism but do not by themselves prove the full PAO/GAO ecological phenotype."
    )
  }),
  Ureolysis = local({
    roles <- aaa_gene_roles(diagnostic = c("ureA", "ureB", "ureC"), supporting = c("ureD", "ureE", "ureF", "ureG"))
    list(
      display_name = "Ureolysis", category = "Nitrogen cycle",
      description = "Evaluates the structural urease complex and accessory maturation genes.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("ureolysis", list(Urease = c("ureA", "ureB", "ureC")), list(Maturation = c("ureD", "ureE", "ureF", "ureG")), 1, .67),
      graph_main = "Ureolytic potential", analysis_name = "Ureolysis",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Sulfide_oxidation = local({
    roles <- aaa_gene_roles(diagnostic = c("sqr", "fccA", "fccB"), supporting = c("soxA", "soxB", "soxX", "soxY", "soxZ", "soxC", "soxD"))
    list(
      display_name = "Sulfide and reduced-sulfur oxidation", category = "Sulfur cycle",
      description = "Evaluates sulfide:quinone oxidoreductase, flavocytochrome c sulfide dehydrogenase and Sox-system evidence.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("reduced-sulfur oxidation", list(Initial_oxidation = c("sqr", "fccA", "fccB")), list(Sox_system = c("soxA", "soxB", "soxX", "soxY", "soxZ", "soxC", "soxD")), .34, .34),
      graph_main = "Reduced-sulfur oxidation potential", analysis_name = "Sulfide_oxidation",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Aerobic_CO_oxidation = local({
    roles <- aaa_gene_roles(diagnostic = c("coxL", "coxM", "coxS"))
    list(
      display_name = "Aerobic carbon monoxide oxidation", category = "Carbon fixation and CO/CO2",
      description = "Evaluates the molybdenum-dependent form-I carbon monoxide dehydrogenase complex.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("aerobic CO oxidation", list(CODH = c("coxL", "coxM", "coxS")), list(), 1, .67),
      graph_main = "Aerobic CO-oxidation potential", analysis_name = "Aerobic_CO_oxidation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "coxL-like sequences require careful annotation because related xanthine-oxidase-family proteins can produce false positives."
    )
  }),
  Carboxydotrophic_hydrogenogenesis = local({
    roles <- aaa_gene_roles(diagnostic = c("cooS", "cooF"), supporting = c("echA", "echB", "echC", "echD", "echE", "echF"))
    list(
      display_name = "Carboxydotrophic hydrogenogenesis (CO to H2)", category = "Fermentation and hydrogen",
      description = "Evaluates anaerobic CO dehydrogenase and energy-converting hydrogenase modules associated with CO-dependent H2 production.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("carboxydotrophic hydrogenogenesis", list(CO_oxidation = c("cooS", "cooF")), list(Ech_hydrogenase = c("echA", "echB", "echC", "echD", "echE", "echF")), 1, .5),
      graph_main = "CO-dependent hydrogenogenic potential", analysis_name = "Carboxydotrophic_hydrogenogenesis",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Fermentative_H2_production = local({
    roles <- aaa_gene_roles(diagnostic = c("hydA", "hydE", "hydF", "hydG"), supporting = c("pflB", "porA", "porB"))
    list(
      display_name = "Fermentative hydrogen production", category = "Fermentation and hydrogen",
      description = "Evaluates [FeFe]-hydrogenase and maturation genes, with fermentative electron-donor modules as supporting evidence.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("fermentative hydrogen production", list(FeFe_hydrogenase = c("hydA"), Maturation = c("hydE", "hydF", "hydG")), list(Fermentative_input = c("pflB", "porA", "porB")), 1, .67),
      graph_main = "Fermentative H2-production potential", analysis_name = "Fermentative_H2_production",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Calvin_cycle = local({
    roles <- aaa_gene_roles(diagnostic = c("rbcL", "rbcS", "prkB"))
    list(
      display_name = "CO2 fixation — Calvin-Benson-Bassham cycle", category = "Carbon fixation and CO/CO2",
      description = "Evaluates Rubisco and phosphoribulokinase evidence for the Calvin cycle.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("Calvin-cycle CO2 fixation", list(Rubisco = c("rbcL", "rbcS"), PRK = c("prkB")), list(), 1, .5),
      graph_main = "Calvin-cycle potential", analysis_name = "Calvin_cycle",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Reverse_TCA_cycle = local({
    roles <- aaa_gene_roles(diagnostic = c("aclA", "aclB"), supporting = c("korA", "korB", "porA", "porB"))
    list(
      display_name = "CO2 fixation — reverse TCA cycle", category = "Carbon fixation and CO/CO2",
      description = "Evaluates ATP-citrate lyase with supporting 2-oxoacid:ferredoxin oxidoreductases.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("reverse-TCA CO2 fixation", list(ATP_citrate_lyase = c("aclA", "aclB")), list(Oxidoreductases = c("korA", "korB", "porA", "porB")), 1, .5),
      graph_main = "Reverse-TCA potential", analysis_name = "Reverse_TCA_cycle",
      pathway_selector = aaa_positive_pathway_call
    )
  }),
  Methylotrophic_methanogenesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("mcrA", "mcrB", "mcrG"),
      supporting = c(
        "mtaA", "mtaB", "mtaC",
        "mtbA", "mtmB", "mtmC", "mtbB", "mtbC", "mttB", "mttC"
      )
    )
    list(
      display_name = "Methylotrophic methanogenesis",
      category = "Methane and C1 metabolism",
      description = paste(
        "Evaluates methylotrophic methanogenesis using a supported Mcr",
        "terminal complex together with a coherent methanol or methylamine",
        "substrate-entry module. Partial modules are reported as insufficient",
        "evidence."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      diagnostic_groups = list(
        Terminal_methanogenesis = c("mcrA", "mcrB", "mcrG"),
        Methanol_entry = c("mtaA", "mtaB", "mtaC"),
        Monomethylamine_entry = c("mtmB", "mtmC", "mtbA"),
        Dimethylamine_entry = c("mtbB", "mtbC", "mtbA"),
        Trimethylamine_entry = c("mttB", "mttC", "mtbA")
      ),
      classifier = classify_methylotrophic_methanogenesis,
      graph_main = "Methylotrophic methanogenesis potential",
      analysis_name = "Methylotrophic_methanogenesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Methyltransferases alone do not establish methanogenesis; the",
        "classifier also requires methanogenesis-specific terminal evidence."
      )
    )
  }),
  Proline_biosynthesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("proA", "proB", "proC")
    )
    list(
      display_name = "Proline biosynthesis",
      category = "Osmoprotectors and high value molecules",
      description = paste(
        "Evaluates the canonical glutamate-to-proline biosynthetic route",
        "through ProB, ProA and ProC."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_proline_biosynthesis,
      graph_main = "Proline biosynthetic potential",
      analysis_name = "Proline_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "The result represents genomic biosynthetic potential.",
        "Proline accumulation and osmoprotection require physiological validation."
      )
    )
  }),
  Pipecolate_biosynthesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("lhpI", "dpkA", "lysDH"),
      supporting = c("proC", "p2cr")
    )
    list(
      display_name = "Pipecolate biosynthesis",
      category = "Osmoprotectors and high value molecules",
      description = paste(
        "Evaluates alternative lysine/pipecolate-associated enzyme routes.",
        "Because several biochemical routes exist, the classifier requires",
        "coherent evidence from entry and reduction steps."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_pipecolate_biosynthesis,
      graph_main = "Pipecolate biosynthetic potential",
      analysis_name = "Pipecolate_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Single markers are reported only as partial evidence.",
        "The prediction does not establish pipecolate production."
      )
    )
  }),
  Carotenoid_pigments = local({
    roles <- aaa_gene_roles(
      diagnostic = c(
        "crtE", "crtB", "crtI", "crtP", "zIso", "crtQ", "crtH",
        "crtY", "crtZ", "crtW"
      ),
      supporting = c(
        "al1", "cruA", "cruP", "al2", "crtR", "lut5",
        "zep", "nsy", "crtC", "crtD", "crtF",
        "crtLe", "crtLb", "lut1", "crtO",
        "carT", "carD", "lyeJ", "cruF"
      )
    )
    list(
      display_name = "Carotenoid pigment production",
      category = "Pigments and phototrophy",
      description = paste(
        "Reconstructs precursor-dependent carotenoid branches and reports",
        "specific pigment potentials only when the required upstream and",
        "product-forming enzymes are coherently detected."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_carotenoid_pigments,
      graph_main = "Carotenoid pigment-production potential",
      analysis_name = "Carotenoid_pigments",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Possible products include lycopene, beta-carotene, zeaxanthin,",
        "violaxanthin, neoxanthin, astaxanthin, canthaxanthin,",
        "spirilloxanthin, lutein, neurosporaxanthin and bacterioruberin.",
        "Predictions indicate genomic potential, not pigment concentration."
      )
    )
  }),
  Bacteriorhodopsin = local({
    roles <- aaa_gene_roles(
      diagnostic = c("bop"),
      supporting = c("blh", "bcmo1")
    )
    list(
      display_name = "Bacteriorhodopsin-based phototrophy",
      category = "Pigments and phototrophy",
      description = paste(
        "Evaluates bacteriorhodopsin together with retinal-producing",
        "beta-carotene cleavage enzymes."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_bacteriorhodopsin,
      graph_main = "Bacteriorhodopsin phototrophic potential",
      analysis_name = "Bacteriorhodopsin",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "A high-confidence call requires bop and at least one retinal",
        "biosynthesis marker. Protein expression and phototrophic activity",
        "are not inferred."
      )
    )
  }),
  Syngas_to_ethanol = local({
    roles <- aaa_gene_roles(
      diagnostic = c(
        "cooS", "acsA", "acsB", "acsC", "acsD", "acsE",
        "adhE", "aor", "adh"
      ),
      supporting = c("cooF", "fhs", "folD", "metF", "fdh")
    )
    list(
      display_name = "Syngas valorisation to ethanol",
      category = "Carbon fixation and CO/CO2",
      description = paste(
        "Evaluates a CO/C1-assimilating Wood-Ljungdahl backbone together",
        "with an ethanol-forming branch. This complements, but does not",
        "duplicate, the existing CO-to-H2 and homoacetogenesis functions."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_syngas_to_ethanol,
      graph_main = "Syngas-to-ethanol valorisation potential",
      analysis_name = "Syngas_to_ethanol",
      pathway_selector = function(x) {
        !grepl("^(No detected|Unknown|C1 assimilation evidence)", x,
          ignore.case = TRUE
        )
      },
      evidence_note = paste(
        "Alcohol dehydrogenases are widespread, so they are interpreted only",
        "when a coherent C1-assimilation backbone is also present.",
        "The result does not predict ethanol yield."
      )
    )
  }),
  Reductive_glycine_pathway = local({
    roles <- aaa_gene_roles(
      diagnostic = c(
        "fdh", "fhs", "folD",
        "gcvP", "gcvT", "gcvH", "lpd", "glyA"
      ),
      supporting = c("sdaA", "sdaB", "tdcG")
    )
    list(
      display_name = "Reductive glycine pathway",
      category = "Carbon fixation and CO/CO2",
      description = paste(
        "Evaluates formate/CO2 assimilation through the tetrahydrofolate",
        "module, reverse glycine-cleavage system, serine formation and",
        "conversion toward pyruvate."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_reductive_glycine_pathway,
      graph_main = "Reductive glycine pathway potential",
      analysis_name = "Reductive_glycine_pathway",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "The glycine-cleavage system is reversible and participates in other",
        "metabolisms. A complete call therefore requires multiple coherent",
        "modules rather than isolated gcv genes."
      )
    )
  }),
  RuMP_formaldehyde_assimilation = local({
    roles <- aaa_gene_roles(
      diagnostic = c("hps", "phi"),
      supporting = c("fae")
    )
    list(
      display_name = "RuMP formaldehyde assimilation",
      category = "Methane and C1 metabolism",
      description = paste(
        "Evaluates the defining hexulose-phosphate synthase and",
        "phosphohexuloisomerase reactions of the ribulose-monophosphate",
        "formaldehyde-assimilation pathway."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_rump_formaldehyde_assimilation,
      graph_main = "RuMP formaldehyde-assimilation potential",
      analysis_name = "RuMP_formaldehyde_assimilation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "This function evaluates formaldehyde assimilation and does not by",
        "itself establish methane oxidation. Methanotrophy remains assessed",
        "by the existing methane-monooxygenase classifier."
      )
    )
  }),
  Dark_fermentation_FHL = local({
    roles <- aaa_gene_roles(
      diagnostic = c("fdhF", "hycE", "hycB", "hycF", "hycG"),
      supporting = c("hycC", "hycD", "hycI", "fhlA")
    )
    list(
      display_name = "Dark fermentation through formate hydrogen lyase",
      category = "Fermentation and hydrogen",
      description = paste(
        "Evaluates conversion of fermentative formate into hydrogen and",
        "carbon dioxide through a coherent FHL/hydrogenase-3 complex."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_dark_fermentation_fhl,
      graph_main = "Formate-hydrogen-lyase potential",
      analysis_name = "Dark_fermentation_FHL",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Hydrogenase or formate-dehydrogenase subunits in isolation are not",
        "treated as a complete fermentative hydrogen-production pathway."
      )
    )
  }),
  Dark_fermentation_PFOR_hydrogenase = local({
    roles <- aaa_gene_roles(
      diagnostic = c(
        "nifJ", "porA", "porB", "porD", "porG",
        "hydA", "hydB", "hydC"
      ),
      supporting = c("fdx")
    )
    list(
      display_name = "PFOR-linked dark-fermentative hydrogen production",
      category = "Fermentation and hydrogen",
      description = paste(
        "Evaluates a pyruvate:ferredoxin oxidoreductase arm coupled to a",
        "ferredoxin-linked [FeFe]-hydrogenase module."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_dark_fermentation_pfor_hydrogenase,
      graph_main = "PFOR-linked dark-fermentative hydrogen potential",
      analysis_name = "Dark_fermentation_PFOR_hydrogenase",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "PFOR and hydrogenases participate in multiple anaerobic processes.",
        "A strong call therefore requires both pathway arms and preferably a",
        "ferredoxin electron carrier."
      )
    )
  }),
  Glycine_betaine_biosynthesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("betA", "betB", "gbsA", "gbsB")
    )
    list(
      display_name = "Glycine-betaine biosynthesis from choline",
      category = "Osmoprotectors and high value molecules",
      description = paste(
        "Evaluates the two-step oxidation of choline through betaine aldehyde",
        "to the compatible solute glycine betaine."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_glycine_betaine_biosynthesis,
      graph_main = "Glycine-betaine biosynthetic potential",
      analysis_name = "Glycine_betaine_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "The classifier accepts BetA/BetB and GbsB/GbsA naming conventions.",
        "Transport and intracellular accumulation are not inferred."
      )
    )
  }),
  Trehalose_biosynthesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("otsA", "otsB", "treY", "treZ"),
      supporting = c("treS")
    )
    list(
      display_name = "Trehalose biosynthesis",
      category = "Osmoprotectors and high value molecules",
      description = paste(
        "Evaluates the independent OtsAB and TreYZ trehalose biosynthetic",
        "routes and reports TreS separately as an interconversion route."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_trehalose_biosynthesis,
      graph_main = "Trehalose biosynthetic potential",
      analysis_name = "Trehalose_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "Trehalose is associated with several stress responses, but pathway",
        "presence does not demonstrate osmotically induced accumulation."
      )
    )
  }),
  Mannosylglycerate_biosynthesis = local({
    roles <- aaa_gene_roles(
      diagnostic = c("mpgS", "mpgP", "mpgSP")
    )
    list(
      display_name = "Mannosylglycerate biosynthesis",
      category = "Osmoprotectors and high value molecules",
      description = paste(
        "Evaluates the two-step MpgS/MpgP route and the bifunctional MpgS/P",
        "route for the compatible solute mannosylglycerate."
      ),
      gene_roles = roles,
      genes = aaa_registry_genes(roles),
      classifier = classify_mannosylglycerate_biosynthesis,
      graph_main = "Mannosylglycerate biosynthetic potential",
      analysis_name = "Mannosylglycerate_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = paste(
        "The prediction indicates biosynthetic capacity and does not establish",
        "stress-dependent synthesis or intracellular concentration."
      )
    )
  }),
  L_lactate_production = local({
    roles <- aaa_gene_roles(diagnostic = c("ldhL"))
    list(
      display_name = "L-lactate production",
      category = "Fermentation and hydrogen",
      description = "Evaluates stereospecific L-lactate dehydrogenase-mediated production potential.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_l_lactate_production,
      graph_main = "L-lactate production potential", analysis_name = "L_lactate_production",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "A stereospecific ldhL annotation is required; generic LDH annotations are not treated as diagnostic."
    )
  }),
  D_lactate_production = local({
    roles <- aaa_gene_roles(diagnostic = c("ldhD"))
    list(
      display_name = "D-lactate production",
      category = "Fermentation and hydrogen",
      description = "Evaluates stereospecific D-lactate dehydrogenase-mediated production potential.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_d_lactate_production,
      graph_main = "D-lactate production potential", analysis_name = "D_lactate_production",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Production is kept separate from respiratory D-lactate oxidation."
    )
  }),
  Respiratory_lactate_utilisation = local({
    roles <- aaa_gene_roles(
      diagnostic = c("lutA", "lutB", "lutC", "lldE", "lldF", "lldG", "lldD", "dld"),
      supporting = c("lutP", "lldP")
    )
    list(
      display_name = "Respiratory lactate utilisation",
      category = "Fermentation and hydrogen",
      description = "Evaluates complete LutABC/LldEFG complexes and stereospecific respiratory lactate dehydrogenases.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_respiratory_lactate_utilisation,
      graph_main = "Respiratory lactate-utilisation potential", analysis_name = "Respiratory_lactate_utilisation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Transporters strengthen the call but do not establish lactate oxidation by themselves."
    )
  }),
  Lactate_to_PHA = local({
    roles <- aaa_gene_roles(
      diagnostic = c("lutA", "lutB", "lutC", "lldE", "lldF", "lldG", "lldD", "dld", "phaA", "phaB", "phaC")
    )
    list(
      display_name = "Lactate-fed PHA biosynthesis",
      category = "Biopolymers",
      description = "Evaluates lactate oxidation coupled to a canonical PhaABC PHA/PHB biosynthetic module.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_lactate_to_pha,
      graph_main = "Lactate-fed PHA potential", analysis_name = "Lactate_to_PHA",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "This predicts use of lactate carbon for PHA metabolism, not direct lactyl-unit incorporation."
    )
  }),
  PLA_engineering_chassis = local({
    roles <- aaa_gene_roles(diagnostic = c("ldhL", "ldhD", "pct", "phaC"))
    list(
      display_name = "Putative PLA engineering chassis",
      category = "Biopolymers",
      description = "Evaluates lactate supply, lactyl-CoA formation and a PHA synthase as prerequisites for engineered lactate-containing polyesters.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_pla_engineering_chassis,
      graph_main = "PLA engineering-chassis potential", analysis_name = "PLA_engineering_chassis",
      pathway_selector = function(x) grepl("Putative", x, fixed = TRUE),
      evidence_note = "Wild-type phaC does not prove PLA polymerisation; enzyme engineering and substrate specificity must be verified experimentally."
    )
  }),
  Methanol_oxidation = local({
    roles <- aaa_gene_roles(
      diagnostic = c("mxaF", "mxaI", "xoxF"),
      supporting = c("mxaG", "mxaJ", "xoxG", "xoxJ")
    )
    list(
      display_name = "Methanol oxidation",
      category = "Methane and C1 metabolism",
      description = "Evaluates PQQ-dependent MxaFI and lanthanide-dependent XoxF methanol dehydrogenase systems.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_methanol_oxidation,
      graph_main = "Methanol-oxidation potential", analysis_name = "Methanol_oxidation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Methanol oxidation is reported separately from methane oxidation and formaldehyde assimilation."
    )
  }),
  Serine_cycle_C1_assimilation = local({
    roles <- aaa_gene_roles(
      diagnostic = c("fae", "mtdA", "mch", "hprA", "sgaA", "mtkA", "mtkB", "mcl"),
      supporting = c("glyA")
    )
    list(
      display_name = "Serine-cycle C1 assimilation",
      category = "Methane and C1 metabolism",
      description = "Evaluates formaldehyde incorporation, serine-to-glycerate conversion and the malyl-CoA cleavage arm.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_serine_cycle_c1_assimilation,
      graph_main = "Serine-cycle C1-assimilation potential", analysis_name = "Serine_cycle_C1_assimilation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Common central enzymes are not sufficient; multiple characteristic modules are required."
    )
  }),
  Glucosylglycerol_biosynthesis = local({
    roles <- aaa_gene_roles(diagnostic = c("ggpS", "ggpP", "ggpSP"))
    list(
      display_name = "Glucosylglycerol biosynthesis",
      category = "Osmoprotectors and high value molecules",
      description = "Evaluates two-step or bifunctional glucosylglycerol biosynthesis.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_glucosylglycerol_biosynthesis,
      graph_main = "Glucosylglycerol biosynthetic potential", analysis_name = "Glucosylglycerol_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Genomic capacity does not demonstrate salt-induced accumulation."
    )
  }),
  NAGGN_biosynthesis = local({
    roles <- aaa_gene_roles(diagnostic = c("asnO", "ngg"))
    list(
      display_name = "NAGGN biosynthesis",
      category = "Osmoprotectors and high value molecules",
      description = "Evaluates the AsnO/Ngg pathway for N-acetylglutaminylglutamine amide biosynthesis.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_naggn_biosynthesis,
      graph_main = "NAGGN biosynthetic potential", analysis_name = "NAGGN_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Both pathway genes are required for a complete call."
    )
  }),
  Sucrose_biosynthesis = local({
    roles <- aaa_gene_roles(diagnostic = c("sps", "spp"))
    list(
      display_name = "Sucrose biosynthesis",
      category = "Osmoprotectors and high value molecules",
      description = "Evaluates the sucrose-phosphate synthase/phosphatase route used by many cyanobacteria as a compatible-solute pathway.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_sucrose_biosynthesis,
      graph_main = "Sucrose biosynthetic potential", analysis_name = "Sucrose_biosynthesis",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "The prediction denotes biosynthetic capacity, not secretion or stress-induced accumulation."
    )
  }),

  # ===========================================================================
  # Biogeochemical cycles.
  # Only cycles whose diagnostic genes are robustly and unambiguously named in
  # RefSeq/GFF annotations are added here. Manganese oxidation, selenium
  # respiration and anaerobic methane oxidation were intentionally deferred:
  # their canonical markers are either poorly conserved or share gene symbols
  # with unrelated pathways, which would produce false positives under the
  # gene-name search Triple_A performs.
  # ===========================================================================

  Phosphonate_degradation = local({
    roles <- aaa_gene_roles(
      diagnostic = c("phnJ"),
      supporting = c("phnG", "phnH", "phnI", "phnK", "phnL", "phnM")
    )
    list(
      display_name = "Phosphonate degradation (C-P lyase)", category = "Phosphorus cycle",
      description = "Evaluates the C-P lyase pathway that cleaves carbon-phosphorus bonds and mobilises phosphate from phosphonates.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("phosphonate degradation", list(CP_lyase = c("phnJ")), list(Phn_cluster = c("phnG", "phnH", "phnI", "phnK", "phnL", "phnM")), 1, .5),
      graph_main = "Phosphonate-degradation potential", analysis_name = "Phosphonate_degradation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "phnJ is the diagnostic radical-SAM C-P lyase subunit; the surrounding phnGHIKLM cluster is supporting evidence."
    )
  }),
  Perchlorate_reduction = local({
    roles <- aaa_gene_roles(diagnostic = c("pcrA", "cld"), supporting = c("pcrB"))
    list(
      display_name = "(Per)chlorate reduction", category = "Halogen cycling",
      description = "Evaluates dissimilatory perchlorate/chlorate reduction using perchlorate reductase and the hallmark chlorite dismutase.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("perchlorate reduction", list(Perchlorate_reductase = c("pcrA"), Chlorite_dismutase = c("cld")), list(Beta_subunit = c("pcrB")), 1, .5),
      graph_main = "(Per)chlorate-reduction potential", analysis_name = "Perchlorate_reduction",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Chlorite dismutase (cld) is the diagnostic hallmark; perchlorate reductase confirms the respiratory route."
    )
  }),
  Organohalide_respiration = local({
    roles <- aaa_gene_roles(diagnostic = c("rdhA"), supporting = c("rdhB"))
    list(
      display_name = "Organohalide respiration (reductive dehalogenation)", category = "Halogen cycling",
      description = "Evaluates reductive dehalogenase evidence for respiratory dehalogenation of chlorinated compounds.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("organohalide respiration", list(Reductive_dehalogenase = c("rdhA")), list(Membrane_anchor = c("rdhB")), 1, .5),
      graph_main = "Organohalide-respiration potential", analysis_name = "Organohalide_respiration",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "rdhA is a diverse gene family; a positive call indicates dehalogenation potential rather than a specific substrate."
    )
  }),
  Dissimilatory_iron_reduction = local({
    roles <- aaa_gene_roles(diagnostic = c("omcA", "omcB", "omcS", "omcZ"))
    list(
      display_name = "Dissimilatory iron reduction", category = "Metals and metalloids",
      description = "Evaluates outer-membrane multiheme cytochromes associated with respiratory Fe(III) reduction (Shewanella and Geobacter lineages).",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_dissimilatory_iron_reduction,
      graph_main = "Dissimilatory iron-reduction potential", analysis_name = "Dissimilatory_iron_reduction",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "Metal-reducing cytochromes are lineage-specific, so the call is capped at 'Probable' and never reports a complete pathway from these markers alone."
    )
  }),
  Iron_oxidation = local({
    roles <- aaa_gene_roles(diagnostic = c("cyc2", "cyc1"))
    list(
      display_name = "Iron (Fe(II)) oxidation", category = "Metals and metalloids",
      description = "Evaluates the outer-membrane cytochrome Cyc2, a diagnostic marker of neutrophilic and acidophilic Fe(II) oxidisers.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = classify_iron_oxidation,
      graph_main = "Iron-oxidation potential", analysis_name = "Iron_oxidation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "cyc2 is the primary diagnostic; the marker indicates oxidation potential, not the specific electron-acceptor coupling."
    )
  }),
  Arsenate_respiration = local({
    roles <- aaa_gene_roles(diagnostic = c("arrA"), supporting = c("arrB"))
    list(
      display_name = "Arsenate respiration", category = "Metals and metalloids",
      description = "Evaluates the respiratory arsenate reductase Arr used for dissimilatory As(V) reduction.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("arsenate respiration", list(Arsenate_reductase = c("arrA")), list(Beta_subunit = c("arrB")), 1, .5),
      graph_main = "Arsenate-respiration potential", analysis_name = "Arsenate_respiration",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "arrA is the catalytic diagnostic subunit; distinguish from the arsenite oxidase aioA."
    )
  }),
  Arsenite_oxidation = local({
    roles <- aaa_gene_roles(diagnostic = c("aioA"), supporting = c("aioB"))
    list(
      display_name = "Arsenite oxidation", category = "Metals and metalloids",
      description = "Evaluates the arsenite oxidase Aio used for As(III) oxidation to less toxic As(V).",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("arsenite oxidation", list(Arsenite_oxidase = c("aioA")), list(Rieske_subunit = c("aioB")), 1, .5),
      graph_main = "Arsenite-oxidation potential", analysis_name = "Arsenite_oxidation",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "aioA is the diagnostic large subunit of arsenite oxidase."
    )
  }),
  Arsenic_detoxification = local({
    roles <- aaa_gene_roles(diagnostic = c("arsC"), supporting = c("arsB", "acr3", "arsA", "arsR"))
    list(
      display_name = "Arsenic detoxification (ars operon)", category = "Metals and metalloids",
      description = "Evaluates the cytoplasmic ars detoxification system (arsenate reductase plus arsenite efflux) conferring arsenic resistance.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("arsenic detoxification", list(Arsenate_reductase = c("arsC")), list(Efflux = c("arsB", "acr3"), ATPase = c("arsA"), Regulator = c("arsR")), 1, .5),
      graph_main = "Arsenic-detoxification potential", analysis_name = "Arsenic_detoxification",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "The ars system is a widespread resistance module and does not by itself indicate dissimilatory arsenic metabolism."
    )
  }),
  Periplasmic_nitrate_reduction = local({
    roles <- aaa_gene_roles(diagnostic = c("napA"), supporting = c("napB", "napC"))
    list(
      display_name = "Periplasmic nitrate reduction (Nap)", category = "Nitrogen cycle",
      description = "Evaluates the periplasmic nitrate reductase Nap, an alternative to the membrane-bound Nar for dissimilatory nitrate reduction.",
      gene_roles = roles, genes = aaa_registry_genes(roles),
      classifier = aaa_make_pathway_classifier("periplasmic nitrate reduction", list(Nitrate_reductase = c("napA")), list(Electron_transfer = c("napB", "napC")), 1, .5),
      graph_main = "Periplasmic nitrate-reduction potential", analysis_name = "Periplasmic_nitrate_reduction",
      pathway_selector = aaa_positive_pathway_call,
      evidence_note = "napA complements the membrane-bound narG evaluated by the denitrification function; it does not by itself imply complete denitrification."
    )
  })
)

aaa_registry <- function() {
  biological_function_registry
}

aaa_registry_functional_definitions <- function(
  ids = names(biological_function_registry)
) {
  missing <- setdiff(ids, names(biological_function_registry))

  if (length(missing) > 0) {
    stop(
      "Unknown biological-function identifiers: ",
      paste(missing, collapse = ", ")
    )
  }

  analyses <- lapply(
    ids,
    function(id) {
      entry <- biological_function_registry[[id]]
      list(
        genes = entry$genes,
        classification_function = entry$classifier,
        graph_main = entry$graph_main,
        graph_note = entry$graph_note %||% NULL,
        analysis_name = entry$analysis_name,
        registry_id = id
      )
    }
  )
  names(analyses) <- ids
  analyses
}

aaa_registry_pathway_definitions <- function(
  ids,
  functional_results = NULL,
  use_analysis_reference = TRUE
) {
  missing <- setdiff(ids, names(biological_function_registry))

  if (length(missing) > 0) {
    stop(
      "Unknown biological-function identifiers: ",
      paste(missing, collapse = ", ")
    )
  }

  pathways <- lapply(
    ids,
    function(id) {
      entry <- biological_function_registry[[id]]
      definition <- list(include = entry$pathway_selector)

      if (use_analysis_reference) {
        definition$analysis <- id
      } else {
        if (is.null(functional_results) ||
          is.null(functional_results[[id]])) {
          stop("Functional result not available for: ", id)
        }
        definition$results <- functional_results[[id]]
      }

      definition
    }
  )

  names(pathways) <- ids
  pathways
}

#' Return one row per gene and curated evidence role
aaa_registry_gene_catalogue <- function(
  ids = names(biological_function_registry)
) {
  missing <- setdiff(ids, names(biological_function_registry))
  if (length(missing) > 0) {
    stop(
      "Unknown biological-function identifiers: ",
      paste(missing, collapse = ", ")
    )
  }

  dplyr::bind_rows(lapply(
    ids,
    function(id) {
      entry <- biological_function_registry[[id]]
      roles <- entry$gene_roles

      dplyr::bind_rows(lapply(
        names(roles),
        function(role) {
          genes <- roles[[role]]
          if (length(genes) == 0) {
            return(NULL)
          }

          tibble::tibble(
            ID = id,
            Function = entry$display_name,
            Category = entry$category,
            Gene = genes,
            Role = role
          )
        }
      ))
    }
  ))
}

#' Return one row per registered biological function
aaa_registry_function_catalogue <- function() {
  dplyr::bind_rows(lapply(
    names(biological_function_registry),
    function(id) {
      entry <- biological_function_registry[[id]]
      roles <- entry$gene_roles

      tibble::tibble(
        ID = id,
        Name = entry$display_name,
        Category = entry$category,
        Description = entry$description,
        Diagnostic_genes = paste(
          roles$diagnostic,
          collapse = ", "
        ),
        Supporting_genes = paste(
          roles$supporting,
          collapse = ", "
        ),
        Accessory_genes = paste(
          roles$accessory,
          collapse = ", "
        ),
        Number_of_genes = length(entry$genes),
        Evidence_note = entry$evidence_note %||% ""
      )
    }
  ))
}

#' Calculate diagnostic/supporting/accessory evidence for classified taxa
aaa_summarize_function_evidence <- function(
  taxa_table,
  function_id
) {
  registry <- aaa_registry()

  if (!function_id %in% names(registry)) {
    stop("Unknown biological function: ", function_id)
  }

  entry <- registry[[function_id]]
  roles <- entry$gene_roles

  missing <- setdiff(entry$genes, names(taxa_table))
  if (length(missing) > 0) {
    stop(
      "The taxon table does not contain registry genes: ",
      paste(missing, collapse = ", ")
    )
  }

  count_role <- function(data, genes) {
    if (length(genes) == 0) {
      return(rep(0L, nrow(data)))
    }
    rowSums(
      as.data.frame(lapply(
        data[genes],
        function(x) x %in% TRUE
      )),
      na.rm = TRUE
    )
  }

  evaluated_role <- function(data, genes) {
    if (length(genes) == 0) {
      return(rep(0L, nrow(data)))
    }
    rowSums(!is.na(data[genes]))
  }

  result <- taxa_table

  for (role in names(roles)) {
    genes <- roles[[role]]
    found_name <- paste0(
      tools::toTitleCase(role),
      "_genes_found"
    )
    evaluated_name <- paste0(
      tools::toTitleCase(role),
      "_genes_evaluated"
    )
    coverage_name <- paste0(
      tools::toTitleCase(role),
      "_coverage"
    )

    found <- count_role(result, genes)
    evaluated <- evaluated_role(result, genes)

    result[[found_name]] <- found
    result[[evaluated_name]] <- evaluated
    result[[coverage_name]] <- ifelse(
      evaluated > 0,
      found / evaluated,
      NA_real_
    )
  }

  # ---------------------------------------------------------------------------
  # Generic module-completeness and confidence tier.
  #
  # Diagnostic genes are treated as the "essential" module for the function;
  # supporting/accessory genes are contextual. Completeness is the KEGG-Module-
  # like fraction of the diagnostic module detected. Confidence additionally
  # accounts for how much of that module could be evaluated at all: a genome
  # whose annotation covers few of the diagnostic genes yields low confidence
  # even if the handful evaluated happen to be present.
  #
  # This is strictly additive. The Potential classification and every existing
  # column are unchanged; only Diagnostic_completeness and Confidence are added.
  # Levels: "High confidence" (>=90% complete), "Medium confidence" (>=50%),
  # "Low confidence" (>0%), "No evidence detected" (0% but adequately assessed),
  # and "Insufficient evidence" (<50% of the diagnostic module could be
  # evaluated, so the call cannot be judged).
  # ---------------------------------------------------------------------------
  core_genes <- if (length(roles$diagnostic) > 0L) roles$diagnostic else entry$genes
  core_total <- length(core_genes)
  core_found <- count_role(result, core_genes)
  core_evaluated <- evaluated_role(result, core_genes)

  completeness <- if (core_total > 0L) {
    core_found / core_total
  } else {
    rep(NA_real_, nrow(result))
  }
  evaluated_fraction <- if (core_total > 0L) {
    core_evaluated / core_total
  } else {
    rep(NA_real_, nrow(result))
  }

  result$Diagnostic_completeness <- round(completeness, 3)

  completeness_safe <- completeness
  completeness_safe[!is.finite(completeness_safe)] <- 0
  insufficient <- !is.finite(evaluated_fraction) |
    core_evaluated == 0L | evaluated_fraction < 0.5

  confidence <- rep("No evidence detected", nrow(result))
  confidence[completeness_safe > 0 & completeness_safe < 0.50] <- "Low confidence"
  confidence[completeness_safe >= 0.50 & completeness_safe < 0.90] <- "Medium confidence"
  confidence[completeness_safe >= 0.90] <- "High confidence"
  confidence[insufficient] <- "Insufficient evidence"
  result$Confidence <- confidence

  if (is.function(entry$evidence_scorer)) {
    scored <- lapply(seq_len(nrow(result)), function(i) {
      marker_values <- vapply(entry$genes, function(gene) {
        value <- result[[gene]][i]
        if (length(value) == 0L || is.na(value)) NA else isTRUE(as.logical(value))
      }, logical(1))
      names(marker_values) <- entry$genes
      entry$evidence_scorer(marker_values)
    })

    result$Evidence_score <- vapply(
      scored,
      function(x) if (is.null(x$score)) NA_real_ else as.numeric(x$score),
      numeric(1)
    )
    result$Evidence_level <- vapply(
      scored,
      function(x) if (is.null(x$level)) "Unknown" else as.character(x$level),
      character(1)
    )
    result$Marker_evaluation_fraction <- vapply(
      scored,
      function(x) {
        if (is.null(x$evaluated_fraction)) NA_real_ else as.numeric(x$evaluated_fraction)
      },
      numeric(1)
    )
  }

  result$Functional_registry_ID <- function_id

  result
}
