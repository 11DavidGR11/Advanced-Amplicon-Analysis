qa_register_test(
  "SCI_001", "regression", "critical",
  "Curated biological classifiers reproduce reference marker patterns",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Extensions", "classification_functions.R"), local = FALSE)
    hydrogen <- c(mcrA=TRUE,mcrB=TRUE,mcrG=TRUE,fwdA=TRUE,fwdB=TRUE,fwdC=TRUE,fwdD=TRUE,
                  fmdA=FALSE,fmdB=FALSE,fmdC=FALSE,fmdD=FALSE,ftr=TRUE,mch=TRUE,mtd=TRUE,mer=TRUE,
                  mtrA=TRUE,mtrB=TRUE,mtrC=TRUE,mtrD=TRUE,mtrE=TRUE,mtrF=TRUE,mtrG=TRUE,mtrH=TRUE,
                  ackA=FALSE,pta=FALSE,acs=FALSE,cdhA=FALSE,cdhB=FALSE,cdhC=FALSE,cdhD=FALSE,cdhE=FALSE)
    acetoclastic <- hydrogen; acetoclastic[] <- FALSE
    acetoclastic[c("mcrA","mcrB","mcrG","ackA","pta","cdhA","cdhB","cdhC","cdhD","cdhE")] <- TRUE
    homo <- c(acsA=TRUE, acsB=TRUE, acsC=TRUE, acsD=TRUE, acsE=TRUE, fhs=TRUE, folD=TRUE, metF=TRUE, fdh=TRUE, pta=TRUE, ackA=TRUE)
    qa_expect_true(grepl("Hydrogenotrophic", classify_hydrogenotrophic_methanogenesis(hydrogen)),
                   "Hydrogenotrophic reference pattern was not recognized")
    qa_expect_true(grepl("Acetoclastic", classify_acetoclastic_methanogenesis(acetoclastic)),
                   "Acetoclastic reference pattern was not recognized")
    qa_expect_true(identical(classify_homoacetogenesis(homo), "Homoacetogenesis potential"),
                   "Homoacetogenesis reference pattern was not recognized")
    homo_score <- score_homoacetogenesis(homo)
    qa_expect_true(is.list(homo_score) && identical(homo_score$level, "High") &&
                     isTRUE(homo_score$score >= 80),
                   "Homoacetogenesis score did not preserve named marker indexing")
  }
)

qa_register_test(
  "SCI_002", "regression", "critical",
  "Functional classifiers receive named logical marker vectors",
  function() {
    src <- paste(readLines(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Core", "aaa_functional_potential.R"), warn=FALSE), collapse="\n")
    qa_expect_true(grepl("names(x) <- genes", src, fixed = TRUE),
                   "Functional classification does not preserve marker names")
    qa_expect_true(grepl("Markers_evaluated", src, fixed=TRUE),
                   "Functional evidence diagnostics are missing")
  }
)

qa_register_test(
  "SCI_003", "regression", "critical",
  "New pigment, syngas and C1 classifiers preserve conservative pathway logic",
  function() {
    source(file.path(
      QA_ROOT, "AAApp", "Common", "Engine", "Extensions",
      "classification_functions.R"
    ), local = FALSE)

    pigment <- c(
      crtE=TRUE, crtB=TRUE, crtI=TRUE, crtY=TRUE,
      crtZ=TRUE, crtW=TRUE
    )
    pigment_call <- classify_carotenoid_pigments(pigment)
    qa_expect_true(
      grepl("astaxanthin", pigment_call, ignore.case=TRUE),
      "A coherent astaxanthin marker hierarchy was not recognized"
    )

    incomplete_pigment <- c(crtZ=TRUE, crtW=TRUE)
    qa_expect_true(
      !grepl("astaxanthin", classify_carotenoid_pigments(incomplete_pigment),
             ignore.case=TRUE),
      "Downstream pigment markers were classified without their precursor route"
    )

    syngas <- c(
      cooS=TRUE, cooF=TRUE, acsA=TRUE, acsB=TRUE,
      acsC=TRUE, acsD=TRUE, acsE=FALSE,
      fhs=TRUE, folD=TRUE, metF=TRUE, fdh=TRUE,
      adhE=TRUE, aor=FALSE, adh=FALSE
    )
    qa_expect_true(
      identical(
        classify_syngas_to_ethanol(syngas),
        "Syngas-to-ethanol valorisation potential"
      ),
      "The coherent syngas-to-ethanol reference pattern was not recognized"
    )

    rgly <- c(
      fdh=TRUE, fhs=TRUE, folD=TRUE,
      gcvP=TRUE, gcvT=TRUE, gcvH=TRUE, lpd=TRUE,
      glyA=TRUE, sdaA=TRUE, sdaB=FALSE, tdcG=FALSE
    )
    qa_expect_true(
      identical(
        classify_reductive_glycine_pathway(rgly),
        "Complete reductive glycine pathway potential"
      ),
      "The reductive glycine reference pattern was not recognized"
    )
  }
)

qa_register_test(
  "SCI_004", "regression", "critical",
  "Methane assimilation, dark fermentation and osmoprotection classifiers require coherent modules",
  function() {
    source(file.path(
      QA_ROOT, "AAApp", "Common", "Engine", "Extensions",
      "classification_functions.R"
    ), local = FALSE)

    qa_expect_true(
      grepl(
        "Complete RuMP",
        classify_rump_formaldehyde_assimilation(c(hps=TRUE, phi=TRUE, fae=TRUE)),
        fixed=TRUE
      ),
      "The defining Hps/Phi RuMP pair was not recognized"
    )

    qa_expect_true(
      !grepl(
        "High-confidence",
        classify_dark_fermentation_fhl(c(
          fdhF=TRUE, hycE=TRUE, hycB=FALSE, hycF=FALSE, hycG=FALSE,
          hycC=FALSE, hycD=FALSE, hycI=FALSE, fhlA=FALSE
        )),
        fixed=TRUE
      ),
      "FHL was overcalled from only FdhF and HycE"
    )

    qa_expect_true(
      identical(
        classify_dark_fermentation_pfor_hydrogenase(c(
          nifJ=TRUE, porA=FALSE, porB=FALSE, porD=FALSE, porG=FALSE,
          hydA=TRUE, hydB=TRUE, hydC=FALSE, fdx=TRUE
        )),
        "High-confidence PFOR-linked dark-fermentative hydrogen potential"
      ),
      "A coherent PFOR/[FeFe]-hydrogenase pattern was not recognized"
    )

    qa_expect_true(
      identical(
        classify_glycine_betaine_biosynthesis(c(
          betA=TRUE, betB=TRUE, gbsA=FALSE, gbsB=FALSE
        )),
        "Complete glycine-betaine biosynthetic potential from choline"
      ),
      "The canonical BetA/BetB route was not recognized"
    )

    qa_expect_true(
      identical(
        classify_trehalose_biosynthesis(c(
          otsA=TRUE, otsB=TRUE, treY=FALSE, treZ=FALSE, treS=FALSE
        )),
        "Complete OtsAB trehalose biosynthetic potential"
      ),
      "The complete OtsAB route was not recognized"
    )
  }
)

qa_register_test(
  "SCI_005", "regression", "critical",
  "Lactate, PLA, C1 and compatible-solute classifiers avoid single-gene overcalling",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Extensions", "classification_functions.R"), local=FALSE)

    qa_expect_true(identical(classify_l_lactate_production(c(ldhL=TRUE)), "L-lactate production potential"),
                   "L-lactate production was not recognized")
    qa_expect_true(identical(classify_d_lactate_production(c(ldhD=TRUE)), "D-lactate production potential"),
                   "D-lactate production was not recognized")

    lut <- c(lutA=TRUE,lutB=TRUE,lutC=TRUE,lutP=TRUE,lldE=FALSE,lldF=FALSE,lldG=FALSE,lldP=FALSE,lldD=FALSE,dld=FALSE)
    qa_expect_true(grepl("High-confidence", classify_respiratory_lactate_utilisation(lut), fixed=TRUE),
                   "Complete LutABC plus transporter was not recognized")

    qa_expect_true(!grepl("Putative", classify_pla_engineering_chassis(c(ldhL=TRUE,pct=TRUE,phaC=FALSE,ldhD=FALSE)), fixed=TRUE),
                   "PLA chassis was overcalled without a polymerase")
    qa_expect_true(grepl("Putative", classify_pla_engineering_chassis(c(ldhL=TRUE,pct=TRUE,phaC=TRUE,ldhD=FALSE)), fixed=TRUE),
                   "Coherent PLA engineering chassis was not recognized")

    serine <- c(fae=TRUE,mtdA=TRUE,mch=TRUE,glyA=TRUE,hprA=TRUE,sgaA=TRUE,mtkA=TRUE,mtkB=FALSE,mcl=TRUE)
    qa_expect_true(identical(classify_serine_cycle_c1_assimilation(serine), "Complete serine-cycle C1-assimilation potential"),
                   "Coherent serine-cycle modules were not recognized")

    qa_expect_true(identical(classify_glucosylglycerol_biosynthesis(c(ggpS=TRUE,ggpP=TRUE,ggpSP=FALSE)),
                             "Complete glucosylglycerol biosynthetic potential"),
                   "GgpS/GgpP route was not recognized")
    qa_expect_true(identical(classify_naggn_biosynthesis(c(asnO=TRUE,ngg=TRUE)), "Complete NAGGN biosynthetic potential"),
                   "AsnO/Ngg route was not recognized")
    qa_expect_true(identical(classify_sucrose_biosynthesis(c(sps=TRUE,spp=TRUE)), "Complete SPS/SPP sucrose biosynthetic potential"),
                   "SPS/SPP route was not recognized")
  }
)



qa_register_test(
  "SCI_006", "regression", "critical",
  "Methanogenesis profile applies route completeness and methanogenic priority",
  function() {
    source(file.path(QA_ROOT, "AAApp", "Common", "Engine", "Extensions", "classification_functions.R"), local = FALSE)

    methyl <- c(mcrA=TRUE,mcrB=TRUE,mcrG=TRUE,mtaA=TRUE,mtaB=TRUE,mtaC=TRUE)
    qa_expect_true(
      identical(classify_methylotrophic_methanogenesis(methyl),
                "Methylotrophic methanogenesis potential"),
      "Independent methylotrophic methanogenesis was not recognized"
    )
    qa_expect_true(
      grepl("Methylotrophic methanogenesis", classify_methanogenesis_subtype(methyl), fixed=TRUE),
      "Combined profile did not recognize methylotrophic methanogenesis"
    )

    homo <- c(acsA=TRUE,acsB=TRUE,acsC=TRUE,acsD=TRUE,acsE=FALSE,
              fhs=TRUE,folD=TRUE,metF=FALSE,fdh=TRUE)
    qa_expect_true(
      identical(classify_methanogenesis_subtype(homo),
                "Reductive Wood-Ljungdahl pathway potential (possible homoacetogenesis)"),
      "Complete non-methanogenic WLP was not recognized"
    )

    methanogenic_wlp <- c(homo, mcrA=TRUE, mcrB=TRUE, mcrG=TRUE,
                          ftr=TRUE, mch=TRUE, mtd=TRUE, mer=TRUE,
                          fwdA=TRUE, fwdB=TRUE, fwdC=TRUE,
                          mtrA=TRUE, mtrB=TRUE, mtrC=TRUE, mtrD=TRUE,
                          mtrE=TRUE, mtrF=TRUE)
    priority_call <- classify_methanogenesis_subtype(methanogenic_wlp)
    qa_expect_true(
      grepl("methanogenesis", priority_call, ignore.case=TRUE) &&
        !grepl("homoacetogenesis", priority_call, ignore.case=TRUE),
      "A methanogenic genome was incorrectly classified as homoacetogenic"
    )
    qa_expect_true(
      identical(classify_homoacetogenesis(methanogenic_wlp),
                "No homoacetogenesis assigned: methanogenic pathway detected"),
      "Independent homoacetogenesis classifier did not apply methanogenic exclusion"
    )

    partial <- c(mcrA=TRUE,mcrB=FALSE,mcrG=FALSE,mtaA=TRUE,mtaB=FALSE)
    qa_expect_true(
      identical(classify_methanogenesis_subtype(partial),
                "Insufficient evidence to classify metabolic potential"),
      "Partial pathways were not grouped as insufficient evidence"
    )

    methyl_without_mcr <- c(mcrA=FALSE,mcrB=FALSE,mcrG=FALSE,mtaA=TRUE,mtaB=TRUE,mtaC=TRUE)
    qa_expect_true(
      !grepl("Methylotrophic methanogenesis potential",
             classify_methanogenesis_subtype(methyl_without_mcr), fixed=TRUE),
      "Methyltransferases were overcalled without the Mcr terminal complex"
    )
  }
)
