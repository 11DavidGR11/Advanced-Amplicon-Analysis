# Gene annotation aliases used by aaa_search_gene().
# The original gene symbol is always searched automatically.

gene_aliases <- list(
  ectA = c(
    "diaminobutyrate acetyltransferase",
    "L-2,4-diaminobutyrate acetyltransferase"
  ),
  ectB = c(
    "diaminobutyrate transaminase",
    "L-2,4-diaminobutyrate transaminase"
  ),
  ectC = c("ectoine synthase"),
  ectD = c("ectoine hydroxylase", "ectoine dioxygenase"),
  pmoA = c("particulate methane monooxygenase subunit A"),
  pmoB = c("particulate methane monooxygenase subunit B"),
  pmoC = c("particulate methane monooxygenase subunit C"),
  mmoX = c(
    "methane monooxygenase hydroxylase alpha subunit",
    "soluble methane monooxygenase component A alpha subunit"
  ),
  mmoY = c(
    "methane monooxygenase hydroxylase beta subunit",
    "soluble methane monooxygenase component A beta subunit"
  ),
  mmoZ = c(
    "methane monooxygenase hydroxylase gamma subunit",
    "soluble methane monooxygenase component A gamma subunit"
  ),
  mmoB = c(
    "methane monooxygenase regulatory protein B",
    "soluble methane monooxygenase regulatory protein"
  ),
  mmoC = c(
    "methane monooxygenase reductase",
    "soluble methane monooxygenase reductase"
  ),
  mmoD = c(
    "methane monooxygenase component D",
    "soluble methane monooxygenase protein D"
  ),
  mcrA = c(
    "methyl-coenzyme M reductase subunit alpha",
    "methyl coenzyme M reductase subunit alpha"
  ),
  mcrB = c(
    "methyl-coenzyme M reductase subunit beta",
    "methyl coenzyme M reductase subunit beta"
  ),
  mcrG = c(
    "methyl-coenzyme M reductase subunit gamma",
    "methyl coenzyme M reductase subunit gamma"
  ),
  amoA = c("ammonia monooxygenase subunit A"),
  amoB = c("ammonia monooxygenase subunit B"),
  amoC = c("ammonia monooxygenase subunit C"),
  hao = c(
    "hydroxylamine dehydrogenase",
    "hydroxylamine oxidoreductase"
  ),
  nxrA = c(
    "nitrite oxidoreductase subunit alpha",
    "nitrite dehydrogenase subunit alpha"
  ),
  nxrB = c(
    "nitrite oxidoreductase subunit beta",
    "nitrite dehydrogenase subunit beta"
  ),
  narG = c(
    "respiratory nitrate reductase alpha subunit",
    "nitrate reductase subunit alpha"
  ),
  nirK = c(
    "copper-containing nitrite reductase",
    "copper nitrite reductase"
  ),
  nirS = c(
    "cytochrome cd1 nitrite reductase",
    "cytochrome cd1-dependent nitrite reductase"
  ),
  norB = c(
    "nitric oxide reductase subunit B",
    "quinol-dependent nitric oxide reductase subunit B"
  ),
  nosZ = c("nitrous oxide reductase"),
  sat = c("sulfate adenylyltransferase", "ATP sulfurylase"),
  aprA = c(
    "adenylylsulfate reductase subunit alpha",
    "APS reductase subunit alpha"
  ),
  aprB = c(
    "adenylylsulfate reductase subunit beta",
    "APS reductase subunit beta"
  ),
  dsrA = c("dissimilatory sulfite reductase subunit alpha"),
  dsrB = c("dissimilatory sulfite reductase subunit beta"),
  nifH = c("nitrogenase iron protein", "nitrogenase reductase"),
  nifD = c(
    "nitrogenase molybdenum-iron protein alpha chain",
    "nitrogenase component I alpha chain"
  ),
  nifK = c(
    "nitrogenase molybdenum-iron protein beta chain",
    "nitrogenase component I beta chain"
  ),
  hzsA = c("hydrazine synthase subunit A"),
  hzsB = c("hydrazine synthase subunit B"),
  hzsC = c("hydrazine synthase subunit C"),
  hdh = c("hydrazine dehydrogenase", "hydrazine oxidoreductase"),
  acsA = c(
    "carbon monoxide dehydrogenase catalytic subunit",
    "anaerobic carbon monoxide dehydrogenase catalytic subunit"
  ),
  acsB = c(
    "acetyl-CoA synthase",
    "carbon monoxide dehydrogenase acetyl-CoA synthase subunit"
  ),
  fhs = c(
    "formate-tetrahydrofolate ligase",
    "10-formyltetrahydrofolate synthetase"
  ),
  fdh = c(
    "formate dehydrogenase",
    "formate dehydrogenase subunit alpha",
    "NAD-dependent formate dehydrogenase",
    "formate dehydrogenase major subunit"
  ),
  acsC = c(
    "corrinoid iron-sulfur protein subunit AcsC",
    "acetyl-CoA synthase corrinoid iron-sulfur protein subunit",
    "carbon monoxide dehydrogenase/acetyl-CoA synthase complex corrinoid protein"
  ),
  acsD = c(
    "corrinoid iron-sulfur protein subunit AcsD",
    "acetyl-CoA synthase corrinoid iron-sulfur protein subunit",
    "carbon monoxide dehydrogenase/acetyl-CoA synthase complex corrinoid protein"
  ),
  acsE = c(
    "5-methyltetrahydrofolate--corrinoid iron-sulfur protein methyltransferase",
    "methyltransferase AcsE",
    "methyltetrahydrofolate:corrinoid iron-sulfur protein methyltransferase"
  ),
  folD = c(
    "methylenetetrahydrofolate dehydrogenase/cyclohydrolase",
    "bifunctional methylenetetrahydrofolate dehydrogenase cyclohydrolase"
  ),
  metF = c(
    "5,10-methylenetetrahydrofolate reductase",
    "methylenetetrahydrofolate reductase"
  ),

  # Hydrogenotrophic methanogenesis: CO2 reduction and methyl transfer
  fwdA = c(
    "formylmethanofuran dehydrogenase subunit A",
    "tungsten-containing formylmethanofuran dehydrogenase subunit A"
  ),
  fwdB = c(
    "formylmethanofuran dehydrogenase subunit B",
    "tungsten-containing formylmethanofuran dehydrogenase subunit B"
  ),
  fwdC = c(
    "formylmethanofuran dehydrogenase subunit C",
    "tungsten-containing formylmethanofuran dehydrogenase subunit C"
  ),
  fwdD = c(
    "formylmethanofuran dehydrogenase subunit D",
    "tungsten-containing formylmethanofuran dehydrogenase subunit D"
  ),
  fmdA = c(
    "formylmethanofuran dehydrogenase subunit A",
    "molybdenum-containing formylmethanofuran dehydrogenase subunit A"
  ),
  fmdB = c(
    "formylmethanofuran dehydrogenase subunit B",
    "molybdenum-containing formylmethanofuran dehydrogenase subunit B"
  ),
  fmdC = c(
    "formylmethanofuran dehydrogenase subunit C",
    "molybdenum-containing formylmethanofuran dehydrogenase subunit C"
  ),
  fmdD = c(
    "formylmethanofuran dehydrogenase subunit D",
    "molybdenum-containing formylmethanofuran dehydrogenase subunit D"
  ),
  ftr = c(
    "formylmethanofuran--tetrahydromethanopterin N-formyltransferase",
    "formylmethanofuran tetrahydromethanopterin formyltransferase"
  ),
  mch = c(
    "methenyltetrahydromethanopterin cyclohydrolase",
    "methenyl-H4MPT cyclohydrolase"
  ),
  mtd = c(
    "methylenetetrahydromethanopterin dehydrogenase",
    "methylene-H4MPT dehydrogenase"
  ),
  mer = c(
    "methylenetetrahydromethanopterin reductase",
    "methylene-H4MPT reductase"
  ),
  mtrA = c("tetrahydromethanopterin S-methyltransferase subunit A"),
  mtrB = c("tetrahydromethanopterin S-methyltransferase subunit B"),
  mtrC = c("tetrahydromethanopterin S-methyltransferase subunit C"),
  mtrD = c("tetrahydromethanopterin S-methyltransferase subunit D"),
  mtrE = c("tetrahydromethanopterin S-methyltransferase subunit E"),
  mtrF = c("tetrahydromethanopterin S-methyltransferase subunit F"),
  mtrG = c("tetrahydromethanopterin S-methyltransferase subunit G"),
  mtrH = c("tetrahydromethanopterin S-methyltransferase subunit H"),

  # Acetoclastic methanogenesis: acetate activation and CODH/ACS complex
  ackA = c(
    "acetate kinase",
    "acetate kinase AckA"
  ),
  pta = c(
    "phosphate acetyltransferase",
    "phosphotransacetylase"
  ),
  acs = c(
    "acetyl-CoA synthetase",
    "acetate--CoA ligase",
    "AMP-forming acetyl-CoA synthetase",
    "acetyl-coenzyme A synthetase"
  ),
  cdhA = c(
    "acetyl-CoA decarbonylase/synthase complex subunit alpha",
    "carbon monoxide dehydrogenase subunit alpha",
    "CODH/ACS complex subunit alpha"
  ),
  cdhB = c(
    "acetyl-CoA decarbonylase/synthase complex subunit beta",
    "carbon monoxide dehydrogenase subunit beta",
    "CODH/ACS complex subunit beta"
  ),
  cdhC = c(
    "acetyl-CoA decarbonylase/synthase complex subunit gamma",
    "CODH/ACS complex subunit gamma"
  ),
  cdhD = c(
    "acetyl-CoA decarbonylase/synthase complex subunit delta",
    "CODH/ACS complex subunit delta"
  ),
  cdhE = c(
    "acetyl-CoA decarbonylase/synthase complex subunit epsilon",
    "CODH/ACS complex subunit epsilon"
  ),
  # Triple_A wastewater, C1 and hydrogen functions
  nrfA = c("cytochrome c nitrite reductase catalytic subunit", "ammonia-forming nitrite reductase NrfA"),
  nrfH = c("cytochrome c nitrite reductase electron transfer subunit NrfH"),
  napA = c("periplasmic nitrate reductase catalytic subunit NapA"),
  nosD = c("nitrous oxide reductase maturation protein NosD"),
  nosF = c("nitrous oxide reductase maturation protein NosF"),
  nosY = c("nitrous oxide reductase maturation protein NosY"),
  nosL = c("nitrous oxide reductase accessory protein NosL"),
  ppk1 = c("polyphosphate kinase 1", "polyphosphate kinase PPK1"),
  ppk2 = c("polyphosphate kinase 2", "polyphosphate kinase PPK2"),
  ppx = c("exopolyphosphatase", "exopolyphosphatase PPX"),
  phoU = c("phosphate signaling protein PhoU", "phosphate transport system regulatory protein PhoU"),
  ureA = c("urease subunit gamma"), ureB = c("urease subunit beta"), ureC = c("urease subunit alpha"),
  ureD = c("urease accessory protein UreD"), ureE = c("urease accessory protein UreE"),
  ureF = c("urease accessory protein UreF"), ureG = c("urease accessory protein UreG"),
  sqr = c("sulfide:quinone oxidoreductase", "sulfide quinone oxidoreductase"),
  fccA = c("flavocytochrome c sulfide dehydrogenase flavoprotein subunit"),
  fccB = c("flavocytochrome c sulfide dehydrogenase cytochrome subunit"),
  soxA = c("L-cysteine S-thiosulfotransferase SoxA", "SoxAX cytochrome complex subunit A"),
  soxB = c("SoxB thiosulfohydrolase"), soxX = c("SoxAX cytochrome complex subunit X"),
  soxY = c("sulfur oxidation protein SoxY"), soxZ = c("sulfur oxidation protein SoxZ"),
  soxC = c("sulfane dehydrogenase subunit SoxC"), soxD = c("sulfane dehydrogenase subunit SoxD"),
  coxL = c("aerobic carbon monoxide dehydrogenase large subunit", "carbon monoxide dehydrogenase CoxL"),
  coxM = c("aerobic carbon monoxide dehydrogenase medium subunit"),
  coxS = c("aerobic carbon monoxide dehydrogenase small subunit"),
  cooS = c("anaerobic carbon monoxide dehydrogenase catalytic subunit CooS"),
  cooF = c("carbon monoxide dehydrogenase iron-sulfur protein CooF"),
  echA = c("energy-converting hydrogenase subunit A"), echB = c("energy-converting hydrogenase subunit B"),
  echC = c("energy-converting hydrogenase subunit C"), echD = c("energy-converting hydrogenase subunit D"),
  echE = c("energy-converting hydrogenase subunit E"), echF = c("energy-converting hydrogenase subunit F"),
  hydA = c("[FeFe] hydrogenase catalytic subunit", "[FeFe]-hydrogenase catalytic subunit", "iron-only hydrogenase", "ferredoxin hydrogenase HydA"),
  hydE = c("[FeFe] hydrogenase maturase HydE"), hydF = c("[FeFe] hydrogenase maturase HydF"),
  hydG = c("[FeFe] hydrogenase maturase HydG"),
  pflB = c("pyruvate formate-lyase"), porA = c("pyruvate:ferredoxin oxidoreductase subunit alpha", "pyruvate:ferredoxin oxidoreductase alpha subunit"),
  porB = c("pyruvate:ferredoxin oxidoreductase subunit beta", "pyruvate:ferredoxin oxidoreductase beta subunit"),
  rbcL = c("ribulose-bisphosphate carboxylase large subunit", "Rubisco large subunit"),
  rbcS = c("ribulose-bisphosphate carboxylase small subunit", "Rubisco small subunit"),
  prkB = c("phosphoribulokinase"), aclA = c("ATP citrate lyase alpha subunit"),
  aclB = c("ATP citrate lyase beta subunit"), korA = c("2-oxoglutarate:ferredoxin oxidoreductase subunit alpha"),
  korB = c("2-oxoglutarate:ferredoxin oxidoreductase subunit beta"),
  mtaA = c("methylcobamide:coenzyme M methyltransferase MtaA"),
  mtaB = c("methanol:corrinoid methyltransferase MtaB"), mtaC = c("methanol corrinoid protein MtaC"),
  mtbA = c("methylcobamide:coenzyme M methyltransferase MtbA"),
  mtmB = c("monomethylamine methyltransferase MtmB"), mtmC = c("monomethylamine corrinoid protein MtmC"),
  mtbB = c("dimethylamine methyltransferase MtbB"), mtbC = c("dimethylamine corrinoid protein MtbC"),
  mttB = c("trimethylamine methyltransferase MttB"), mttC = c("trimethylamine corrinoid protein MttC"),
  # Additional non-duplicated compatible-solute, pigment and C1 functions
  proA = c(
    "gamma-glutamyl phosphate reductase",
    "glutamate-5-semialdehyde dehydrogenase"
  ),
  proB = c(
    "glutamate 5-kinase",
    "glutamate kinase"
  ),
  proC = c(
    "pyrroline-5-carboxylate reductase",
    "delta-1-pyrroline-5-carboxylate reductase"
  ),
  lhpI = c(
    "L-pipecolate dehydrogenase",
    "pipecolate oxidoreductase"
  ),
  dpkA = c(
    "delta-1-piperideine-2-carboxylate reductase",
    "D-lysine 5,6-aminomutase-associated pipecolate reductase",
    "L-pipecolate reductase"
  ),
  lysDH = c(
    "lysine dehydrogenase",
    "L-lysine dehydrogenase"
  ),
  p2cr = c(
    "pyrroline-2-carboxylate reductase",
    "delta-1-pyrroline-2-carboxylate reductase"
  ),
  crtE = c(
    "geranylgeranyl diphosphate synthase",
    "geranylgeranyl pyrophosphate synthase"
  ),
  crtB = c(
    "phytoene synthase",
    "15-cis-phytoene synthase"
  ),
  crtP = c("phytoene desaturase CrtP"),
  zIso = c(
    "15-cis-zeta-carotene isomerase",
    "zeta-carotene isomerase"
  ),
  crtQ = c("zeta-carotene desaturase"),
  crtH = c("prolycopene isomerase", "carotenoid isomerase CrtH"),
  crtI = c("phytoene desaturase CrtI"),
  al1 = c("phytoene desaturase AL-1", "carotenoid desaturase AL-1"),
  crtY = c(
    "lycopene beta-cyclase",
    "lycopene cyclase CrtY"
  ),
  cruA = c("lycopene cyclase CruA"),
  cruP = c("lycopene cyclase CruP"),
  al2 = c("lycopene cyclase AL-2"),
  crtR = c("beta-carotene hydroxylase CrtR"),
  crtZ = c(
    "beta-carotene 3-hydroxylase",
    "beta-carotene hydroxylase CrtZ"
  ),
  lut5 = c("carotenoid beta-ring hydroxylase LUT5"),
  zep = c(
    "zeaxanthin epoxidase",
    "zeaxanthin epoxidase ABA1"
  ),
  nsy = c("neoxanthin synthase"),
  crtW = c(
    "beta-carotene ketolase",
    "carotenoid ketolase CrtW"
  ),
  crtC = c("hydroxyneurosporene synthase", "carotenoid hydratase CrtC"),
  crtD = c(
    "methoxyneurosporene dehydrogenase",
    "carotenoid 3,4-desaturase CrtD"
  ),
  crtF = c("hydroxyneurosporene-O-methyltransferase", "carotenoid methyltransferase CrtF"),
  crtLe = c("lycopene epsilon-cyclase", "lycopene epsilon cyclase"),
  crtLb = c("lycopene beta-cyclase", "lycopene beta cyclase"),
  lut1 = c("carotenoid epsilon-ring hydroxylase LUT1"),
  crtO = c("beta-carotene ketolase CrtO", "beta-carotene 4-ketolase"),
  carT = c("carotenoid cleavage dioxygenase CarT"),
  carD = c("carotenoid aldehyde dehydrogenase CarD"),
  lyeJ = c("lycopene elongase LyeJ", "bacterioruberin synthase LyeJ"),
  cruF = c("bisanhydrobacterioruberin hydratase CruF"),
  bcmo1 = c(
    "beta-carotene 15,15'-monooxygenase",
    "beta-carotene 15,15'-dioxygenase"
  ),
  blh = c(
    "beta-carotene 15,15'-dioxygenase Blh",
    "beta-carotene cleavage enzyme Blh"
  ),
  bop = c(
    "bacteriorhodopsin",
    "bacterio-opsin"
  ),
  adhE = c(
    "bifunctional acetaldehyde-CoA/alcohol dehydrogenase",
    "aldehyde-alcohol dehydrogenase AdhE"
  ),
  aor = c(
    "aldehyde:ferredoxin oxidoreductase",
    "aldehyde ferredoxin oxidoreductase"
  ),
  adh = c(
    "alcohol dehydrogenase",
    "NAD-dependent alcohol dehydrogenase"
  ),
  gcvP = c(
    "glycine dehydrogenase decarboxylating",
    "glycine cleavage system P protein"
  ),
  gcvT = c(
    "aminomethyltransferase",
    "glycine cleavage system T protein"
  ),
  gcvH = c(
    "glycine cleavage system H protein",
    "glycine cleavage protein H"
  ),
  lpd = c(
    "dihydrolipoyl dehydrogenase",
    "dihydrolipoamide dehydrogenase"
  ),
  glyA = c(
    "serine hydroxymethyltransferase",
    "glycine hydroxymethyltransferase"
  ),
  sdaA = c("L-serine dehydratase SdaA", "serine deaminase SdaA"),
  sdaB = c("L-serine dehydratase SdaB", "serine deaminase SdaB"),
  tdcG = c("L-serine dehydratase TdcG", "serine deaminase TdcG"),
  # Methane-derived formaldehyde assimilation
  hps = c(
    "3-hexulose-6-phosphate synthase",
    "3-hexulose 6-phosphate synthase",
    "hexulose phosphate synthase"
  ),
  phi = c(
    "6-phospho-3-hexuloisomerase",
    "6-phospho-3-hexulose isomerase",
    "phosphohexuloisomerase",
    "hexulose-6-phosphate isomerase"
  ),
  fae = c(
    "formaldehyde-activating enzyme",
    "formaldehyde activating enzyme"
  ),

  # Formate-hydrogen-lyase dark fermentation
  fdhF = c(
    "formate dehydrogenase-H",
    "formate dehydrogenase H",
    "formate dehydrogenase-H alpha subunit"
  ),
  hycE = c(
    "hydrogenase-3 large subunit",
    "formate hydrogenlyase subunit 5",
    "hydrogenase 3 catalytic subunit"
  ),
  hycB = c(
    "formate hydrogenlyase subunit 2",
    "hydrogenase-3 iron-sulfur subunit HycB"
  ),
  hycF = c(
    "formate hydrogenlyase subunit 6",
    "hydrogenase-3 iron-sulfur subunit HycF"
  ),
  hycG = c(
    "formate hydrogenlyase subunit 7",
    "hydrogenase-3 iron-sulfur subunit HycG"
  ),
  hycC = c("formate hydrogenlyase subunit 3", "hydrogenase-3 subunit HycC"),
  hycD = c("formate hydrogenlyase subunit 4", "hydrogenase-3 subunit HycD"),
  hycI = c(
    "hydrogenase-3 maturation protease",
    "hydrogenase 3 specific protease HycI"
  ),
  fhlA = c(
    "formate hydrogenlyase transcriptional activator",
    "formate hydrogenlyase activator FhlA"
  ),

  # PFOR-linked dark fermentative hydrogen
  nifJ = c(
    "pyruvate:ferredoxin oxidoreductase",
    "pyruvate ferredoxin oxidoreductase NifJ"
  ),
  porD = c("pyruvate:ferredoxin oxidoreductase delta subunit"),
  porG = c("pyruvate:ferredoxin oxidoreductase gamma subunit"),
  hydB = c("[FeFe]-hydrogenase subunit HydB"),
  hydC = c("[FeFe]-hydrogenase subunit HydC"),
  fdx = c("ferredoxin", "2Fe-2S ferredoxin", "4Fe-4S ferredoxin"),

  # Additional compatible-solute routes
  betA = c(
    "choline dehydrogenase",
    "oxygen-dependent choline dehydrogenase"
  ),
  betB = c(
    "betaine-aldehyde dehydrogenase",
    "glycine betaine aldehyde dehydrogenase"
  ),
  gbsA = c(
    "glycine betaine aldehyde dehydrogenase GbsA",
    "betaine-aldehyde dehydrogenase GbsA"
  ),
  gbsB = c(
    "choline dehydrogenase GbsB",
    "type III alcohol dehydrogenase GbsB"
  ),
  otsA = c(
    "trehalose-6-phosphate synthase",
    "alpha,alpha-trehalose-phosphate synthase"
  ),
  otsB = c(
    "trehalose-6-phosphate phosphatase",
    "trehalose-phosphatase"
  ),
  treY = c(
    "maltooligosyltrehalose synthase",
    "maltooligosyl trehalose synthase"
  ),
  treZ = c(
    "maltooligosyltrehalose trehalohydrolase",
    "maltooligosyl trehalose hydrolase"
  ),
  treS = c(
    "maltose alpha-D-glucosyltransferase",
    "trehalose synthase TreS"
  ),
  mpgS = c(
    "mannosyl-3-phosphoglycerate synthase",
    "mannosylglycerate phosphate synthase"
  ),
  mpgP = c(
    "mannosyl-3-phosphoglycerate phosphatase",
    "mannosylglycerate phosphate phosphatase"
  ),
  mpgSP = c(
    "bifunctional mannosyl-3-phosphoglycerate synthase/phosphatase",
    "bifunctional mannosylglycerate synthase phosphatase"
  ),
  # Lactate production and utilisation
  ldhL = c("L-lactate dehydrogenase", "L-lactate dehydrogenase LdhL", "L-specific lactate dehydrogenase"),
  ldhD = c("D-lactate dehydrogenase", "D-lactate dehydrogenase LdhD", "D-specific lactate dehydrogenase"),
  lutA = c("lactate utilization protein LutA", "L-lactate utilization protein A"),
  lutB = c("lactate utilization protein LutB", "L-lactate utilization protein B"),
  lutC = c("lactate utilization protein LutC", "L-lactate utilization protein C"),
  lutP = c("lactate permease LutP", "L-lactate transporter LutP"),
  lldE = c("L-lactate dehydrogenase subunit LldE", "lactate utilization protein LldE"),
  lldF = c("L-lactate dehydrogenase subunit LldF", "lactate utilization protein LldF"),
  lldG = c("L-lactate dehydrogenase subunit LldG", "lactate utilization protein LldG"),
  lldP = c("L-lactate permease", "lactate transporter LldP"),
  lldD = c("L-lactate dehydrogenase LldD", "quinone-dependent L-lactate dehydrogenase"),
  dld = c("D-lactate dehydrogenase Dld", "quinone-dependent D-lactate dehydrogenase"),

  # PHA and engineered lactate-containing polyesters
  phaA = c("acetyl-CoA acetyltransferase PhaA", "beta-ketothiolase PhaA"),
  phaB = c("acetoacetyl-CoA reductase PhaB", "3-ketoacyl-CoA reductase PhaB"),
  phaC = c("polyhydroxyalkanoate synthase", "PHA synthase PhaC", "poly-beta-hydroxybutyrate polymerase"),
  pct = c("propionate CoA-transferase", "propionyl-CoA transferase", "propionate CoA transferase Pct"),

  # Methanol oxidation and serine cycle
  mxaF = c("methanol dehydrogenase large subunit MxaF", "PQQ-dependent methanol dehydrogenase alpha subunit"),
  mxaI = c("methanol dehydrogenase small subunit MxaI", "PQQ-dependent methanol dehydrogenase beta subunit"),
  mxaG = c("methanol dehydrogenase cytochrome c subunit MxaG", "cytochrome c-L MxaG"),
  mxaJ = c("methanol dehydrogenase accessory protein MxaJ"),
  xoxF = c("lanthanide-dependent methanol dehydrogenase XoxF", "PQQ-dependent alcohol dehydrogenase XoxF"),
  xoxG = c("XoxF-associated cytochrome c XoxG"),
  xoxJ = c("XoxF-associated periplasmic protein XoxJ"),
  mtdA = c("methylene-tetrahydromethanopterin dehydrogenase MtdA", "NAD-dependent methylene-H4MPT dehydrogenase"),
  hprA = c("hydroxypyruvate reductase HprA", "glycerate dehydrogenase HprA"),
  sgaA = c("serine-glyoxylate aminotransferase SgaA", "serine glyoxylate transaminase"),
  mtkA = c("malate-CoA ligase alpha subunit MtkA"),
  mtkB = c("malate-CoA ligase beta subunit MtkB"),
  mcl = c("malyl-CoA lyase", "malyl-CoA/beta-methylmalyl-CoA lyase"),

  # Additional compatible solutes
  ggpS = c("glucosylglycerol-phosphate synthase", "glucosyl-glycerol-phosphate synthase GgpS"),
  ggpP = c("glucosylglycerol-phosphate phosphatase", "glucosyl-glycerol-phosphate phosphatase GgpP"),
  ggpSP = c("bifunctional glucosylglycerol-phosphate synthase/phosphatase", "fused GgpS-GgpP protein"),
  asnO = c("NAGGN synthetase AsnO", "glutamine amidotransferase AsnO", "N-acetylglutaminylglutamine amide synthetase"),
  ngg = c("NAGGN acetyltransferase Ngg", "N-acetylglutaminylglutamine amide acetyltransferase"),
  sps = c("sucrose-phosphate synthase", "sucrose-6-phosphate synthase SPS"),
  spp = c("sucrose-phosphate phosphatase", "sucrose-6-phosphate phosphatase SPP"),

  # Biogeochemical cycles
  # Phosphorus — C-P lyase
  phnJ = c("alpha-D-ribose 1-methylphosphonate 5-phosphate C-P lyase", "carbon-phosphorus lyase PhnJ", "C-P lyase subunit PhnJ"),
  phnG = c("carbon-phosphorus lyase complex subunit PhnG"),
  phnH = c("carbon-phosphorus lyase complex subunit PhnH"),
  phnI = c("carbon-phosphorus lyase complex subunit PhnI"),
  phnK = c("phosphonate C-P lyase system protein PhnK", "ABC transporter ATP-binding protein PhnK"),
  phnL = c("phosphonate C-P lyase system protein PhnL", "ABC transporter ATP-binding protein PhnL"),
  phnM = c("alpha-D-ribose 1-methylphosphonate 5-triphosphate diphosphatase", "C-P lyase system protein PhnM"),
  # Chlorine — (per)chlorate reduction and reductive dehalogenation
  pcrA = c("perchlorate reductase subunit alpha", "chlorate reductase subunit alpha PcrA"),
  pcrB = c("perchlorate reductase subunit beta"),
  cld = c("chlorite dismutase", "chlorite O2-lyase Cld"),
  rdhA = c("reductive dehalogenase", "reductive dehalogenase catalytic subunit RdhA"),
  rdhB = c("reductive dehalogenase membrane anchor RdhB", "reductive dehalogenase anchoring protein"),
  # Iron — dissimilatory reduction and Fe(II) oxidation
  omcA = c("outer membrane cytochrome OmcA", "decaheme cytochrome c OmcA"),
  omcB = c("outer membrane cytochrome OmcB", "outer membrane c-type cytochrome OmcB"),
  omcS = c("outer membrane cytochrome OmcS", "multiheme c-type cytochrome OmcS"),
  omcZ = c("outer membrane cytochrome OmcZ", "multiheme c-type cytochrome OmcZ"),
  cyc2 = c("outer membrane cytochrome Cyc2", "iron-oxidizing cytochrome Cyc2", "cytochrome c Cyc2"),
  cyc1 = c("cytochrome c1 Cyc1", "iron oxidation cytochrome Cyc1"),
  # Arsenic — respiration, oxidation and detoxification
  arrA = c("arsenate reductase (respiratory) subunit alpha ArrA", "respiratory arsenate reductase alpha subunit"),
  arrB = c("arsenate reductase (respiratory) subunit beta ArrB", "respiratory arsenate reductase beta subunit"),
  aioA = c("arsenite oxidase large subunit AioA", "arsenite oxidase molybdopterin subunit"),
  aioB = c("arsenite oxidase small subunit AioB", "arsenite oxidase Rieske subunit"),
  arsC = c("arsenate reductase ArsC", "arsenate reductase (glutaredoxin)", "arsenate reductase (thioredoxin)"),
  arsB = c("arsenite efflux transporter ArsB", "arsenical pump membrane protein ArsB"),
  acr3 = c("arsenite efflux transporter Acr3", "arsenical resistance protein Acr3", "arsenite permease Acr3"),
  arsA = c("arsenite-transporting ATPase ArsA", "arsenical pump-driving ATPase"),
  arsR = c("arsenical resistance transcriptional regulator ArsR", "ArsR family transcriptional regulator"),
  # Nitrogen — periplasmic nitrate reduction
  napA = c("periplasmic nitrate reductase NapA", "nitrate reductase catalytic subunit NapA"),
  napB = c("periplasmic nitrate reductase small subunit NapB", "nitrate reductase cytochrome c550-type subunit"),
  napC = c("periplasmic nitrate reductase cytochrome c NapC", "cytochrome c-type protein NapC")
)


# User-maintained aliases created by the Biological Function Builder.
aaa_custom_gene_alias_file <- function(project_root = getOption("triple_a_root", getwd())) {
  file.path(project_root, "Resources", "FunctionalDB", "CustomGeneAliases.json")
}

aaa_normalize_gene_alias_updates <- function(x) {
  if (is.null(x) || !length(x)) {
    return(list())
  }
  if (!is.list(x) || is.null(names(x))) stop("Gene aliases must be a named object.", call. = FALSE)
  output <- list()
  for (gene in names(x)) {
    canonical <- trimws(as.character(gene))
    if (!grepl("^[A-Za-z][A-Za-z0-9_.-]*$", canonical)) stop("Invalid canonical gene name: ", canonical, call. = FALSE)
    aliases <- x[[gene]]
    if (is.list(aliases)) aliases <- unlist(aliases, recursive = TRUE, use.names = FALSE)
    aliases <- unique(trimws(as.character(aliases)))
    aliases <- aliases[nzchar(aliases) & !is.na(aliases) & tolower(aliases) != tolower(canonical)]
    output[[canonical]] <- aliases
  }
  output
}

aaa_validate_gene_alias_dictionary <- function(dictionary) {
  dictionary <- aaa_normalize_gene_alias_updates(dictionary)
  owners <- list()
  for (gene in names(dictionary)) {
    terms <- unique(c(gene, dictionary[[gene]]))
    for (term in terms) {
      key <- tolower(trimws(term))
      previous <- owners[[key]]
      if (!is.null(previous) && !identical(previous, gene)) {
        stop("Alias '", term, "' is assigned to both '", previous, "' and '", gene, "'.", call. = FALSE)
      }
      owners[[key]] <- gene
    }
  }
  invisible(dictionary)
}

aaa_merge_gene_aliases <- function(base, updates) {
  updates <- aaa_normalize_gene_alias_updates(updates)
  aaa_validate_gene_alias_dictionary(updates)
  merged <- base
  for (gene in names(updates)) {
    new_terms <- unique(c(gene, updates[[gene]]))
    other_genes <- setdiff(names(merged), gene)
    for (other in other_genes) {
      occupied <- tolower(unique(c(other, merged[[other]])))
      conflict <- new_terms[tolower(new_terms) %in% occupied]
      if (length(conflict)) {
        stop("Alias '", conflict[[1L]], "' is already assigned to gene '", other, "'.", call. = FALSE)
      }
    }
    merged[[gene]] <- unique(c(merged[[gene]], updates[[gene]]))
  }
  merged
}

aaa_load_custom_gene_aliases <- function(project_root = getOption("triple_a_root", getwd())) {
  file <- aaa_custom_gene_alias_file(project_root)
  if (!file.exists(file)) {
    return(list())
  }
  raw <- jsonlite::read_json(file, simplifyVector = FALSE)
  if (!is.null(raw$genes)) raw <- raw$genes
  aaa_normalize_gene_alias_updates(raw)
}

aaa_save_custom_gene_aliases <- function(updates, project_root = getOption("triple_a_root", getwd())) {
  updates <- aaa_normalize_gene_alias_updates(updates)
  existing <- aaa_load_custom_gene_aliases(project_root)
  combined_custom <- aaa_merge_gene_aliases(existing, updates)
  aaa_merge_gene_aliases(gene_aliases, combined_custom)

  destination <- aaa_custom_gene_alias_file(project_root)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".gene_aliases_", tmpdir = dirname(destination), fileext = ".json")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  jsonlite::write_json(list(genes = combined_custom), temporary, pretty = TRUE, auto_unbox = FALSE)
  if (file.exists(destination) && !file.copy(destination, paste0(destination, ".bak"), overwrite = TRUE)) {
    stop("The previous custom gene dictionary could not be backed up.", call. = FALSE)
  }
  if (!file.rename(temporary, destination) && !file.copy(temporary, destination, overwrite = TRUE)) {
    stop("The custom gene dictionary could not be installed.", call. = FALSE)
  }
  invisible(combined_custom)
}

custom_gene_aliases <- tryCatch(
  aaa_load_custom_gene_aliases(getOption("triple_a_root", getwd())),
  error = function(e) {
    warning("Custom Gene Dictionary could not be loaded: ", conditionMessage(e), call. = FALSE)
    list()
  }
)
gene_aliases <- aaa_merge_gene_aliases(gene_aliases, custom_gene_aliases)
