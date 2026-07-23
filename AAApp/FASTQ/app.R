required_fastq_packages <- c("shiny", "bslib", "dada2", "jsonlite")
# Presence check only: system.file() resolves the install path without loading
# the namespace. dada2 is a heavy Bioconductor package and is used lazily via
# dada2:: once a run starts, so the module UI opens immediately instead of
# blocking on the dada2 load at startup (which could approach the launcher's
# port-wait timeout on slower machines).
missing_fastq_packages <- required_fastq_packages[!nzchar(vapply(required_fastq_packages, function(p) system.file(package = p), character(1)))]
if (length(missing_fastq_packages)) stop("Triple_A FASTQ cannot start. Missing package(s): ", paste(missing_fastq_packages, collapse=", "), call.=FALSE)

library(shiny); library(bslib)
`%||%` <- function(x, y) if (is.null(x)) y else x
app_dir <- normalizePath(getOption("triple_a_fastq_app_dir", getwd()), winslash="/", mustWork=TRUE)
source(file.path(dirname(app_dir), "Common", "paths.R"), local=TRUE)
project_root <- normalizePath(getOption("triple_a_root", aaa_distribution_root(app_dir)), winslash="/", mustWork=TRUE)
fastq_runs_root <- file.path(project_root,"Results","FASTQ"); logs_root <- file.path(project_root,"Results","Logs")
dir.create(fastq_runs_root, recursive=TRUE, showWarnings=FALSE); dir.create(logs_root, recursive=TRUE, showWarnings=FALSE)
source(file.path(project_root,"AAApp","Common","help.R"), local=TRUE)

fmt_time <- function(s) { s <- max(1,round(s)); if(s<60) return(paste0(s," s")); m<-ceiling(s/60); if(m<60) return(paste0(m," min")); paste0(m%/%60," h ",m%%60," min") }
log_line <- function(path, x) cat(sprintf("[%s] %s\n",format(Sys.time(),tz="UTC",usetz=TRUE),x),file=path,append=TRUE)
exe_available <- function(x) nzchar(Sys.which(x))
sample_key <- function(x) {
  # Strip the extension, then ONE read-direction suffix: a separator, an
  # optional "R", the 1/2 direction and an optional lane block. That single
  # pattern already covers _1, .1, -1, _R1 and _R1_001.
  #
  # A second `sub("([._-])[12]$", ...)` used to run here as well, which stripped
  # a further trailing ".1"/".2" — i.e. part of the sample name itself whenever
  # it ended in 1 or 2. Samples such as D90.4.1 and D90.4.2 both collapsed to
  # "D90.4" and pair_files() then aborted with "sample names are ambiguous".
  x <- basename(x); x <- sub("\\.gz$","",x,ignore.case=TRUE); x <- sub("\\.(fastq|fq)$","",x,ignore.case=TRUE)
  x <- sub("([._-])R?[12]([._-]?[0-9]{3})?$","",x,ignore.case=TRUE)
  x
}
pair_files <- function(fwd_names, rev_names) {
  fk <- sample_key(fwd_names); rk <- sample_key(rev_names)
  if (anyDuplicated(fk) || anyDuplicated(rk)) stop("Paired-end sample names are ambiguous after removing R1/R2 suffixes.")
  idx <- match(fk,rk); if(anyNA(idx) || length(fk)!=length(rk)) stop("Forward and reverse FASTQ files could not be paired by sample name.")
  list(samples=make.unique(fk), reverse_order=idx)
}
write_fasta <- function(seqs,path){ con<-file(path,"wt"); on.exit(close(con),add=TRUE); for(i in seq_along(seqs)) writeLines(c(paste0(">ASV",i),seqs[[i]]),con) }
run_command <- function(command,args,log){ log_line(log,paste(command,paste(args,collapse=" "))); status<-system2(command,args,stdout=log,stderr=log); if(status!=0) stop(command," failed with exit status ",status) }

ui <- page_sidebar(
  title = "FASTQ Pipeline (Triple A)",
  theme=bs_theme(version=5,bootswatch="flatly",primary="#5B3F7A"),
  sidebar=sidebar(width=390,
    aaa_help_button(), tags$hr(),
    tags$div(class="alert alert-info small",tags$strong("Scope: "),"technical FASTQ preprocessing and ASV generation only."),
    selectInput("mode","Read layout",c("Single-end"="single","Paired-end"="paired")),
    fileInput("forward","FASTQ / forward (R1) files",multiple=TRUE,accept=c(".fastq",".fq",".gz")),
    conditionalPanel("input.mode == 'paired'",fileInput("reverse","Reverse (R2) files",multiple=TRUE,accept=c(".fastq",".fq",".gz"))),
    fileInput("training","Taxonomy training FASTA",accept=c(".fa",".fasta",".gz")),
    accordion(open=FALSE,
      accordion_panel("Optional external preprocessing",
        checkboxInput("use_fastp","Use fastp for QC/adapters when installed",FALSE),
        checkboxInput("use_cutadapt","Use cutadapt for explicit primer removal when installed",FALSE),
        textInput("primer_f","Forward primer (5'→3')",""),
        textInput("primer_r","Reverse primer (5'→3')",""),
        tags$p(class="small text-muted","External tools are optional. Their availability is checked before the run."),
        tags$p(class="small text-muted",
          "Not installed on Windows? cutadapt: install Python (winget install Python.Python.3.12) then 'pip install cutadapt'. ",
          "fastp has no native Windows build; run it via WSL/Conda (bioconda) or Docker. ",
          "See the contextual FASTQ help for step-by-step instructions.")),
      accordion_panel("DADA2 filtering",
        numericInput("trunc_f","Forward truncation length (0 = none)",240,min=0),
        conditionalPanel("input.mode == 'paired'",numericInput("trunc_r","Reverse truncation length (0 = none)",200,min=0)),
        numericInput("max_ee_f","Forward max expected errors",2,min=0,step=.5),
        conditionalPanel("input.mode == 'paired'",numericInput("max_ee_r","Reverse max expected errors",2,min=0,step=.5)),
        numericInput("min_len","Minimum length",50,min=1),
        numericInput("min_overlap","Minimum paired overlap",12,min=1),
        checkboxInput("pool","Pseudo-pool samples",FALSE))
    ),
    textInput("run_name","Run name","FASTQ_Run"),
    actionButton("run","Process FASTQ",class="btn-primary w-100")
  ),
  navset_card_tab(
    nav_panel("Overview",layout_columns(value_box("Status",textOutput("status")),value_box("External tools",uiOutput("tools")),value_box("Run folder",textOutput("folder")),col_widths=c(5,3,4)),card(card_header("Technical summary"),tableOutput("summary"))),
    nav_panel("Read tracking",card(plotOutput("tracking_plot",height=420),tableOutput("tracking"))),
    nav_panel("ASV table",card(tableOutput("asv"))),
    nav_panel("Taxonomy",card(tableOutput("taxa"))),
    nav_panel("Files and log",layout_columns(card(card_header("Generated files"),uiOutput("files")),card(card_header("Log"),verbatimTextOutput("log")),col_widths=c(5,7)))
  )
)

server <- function(input,output,session){
  aaa_register_context_help(input, output, session, project_root, "FASTQ")
  rv<-reactiveValues(status="Ready",run_dir=NULL,log=NULL,tracking=NULL,asv=NULL,taxa=NULL,summary=NULL)
  output$status<-renderText(rv$status); output$folder<-renderText(if(is.null(rv$run_dir)) "Not created" else rv$run_dir)
  output$tools<-renderUI(tags$span(paste0("fastp: ",if(exe_available("fastp"))"yes" else "no"," · cutadapt: ",if(exe_available("cutadapt"))"yes" else "no")))
  output$summary<-renderTable(rv$summary,striped=TRUE); output$tracking<-renderTable(rv$tracking,striped=TRUE)
  output$tracking_plot<-renderPlot({ req(rv$tracking); x<-rv$tracking; matplot(t(as.matrix(x[-1])),type="b",pch=19,lty=1,xaxt="n",ylab="Reads",xlab="Stage"); axis(1,1:(ncol(x)-1),names(x)[-1]); legend("topright",legend=x$sample,lty=1,pch=19,cex=.7) })
  output$asv<-renderTable({req(rv$asv);head(rv$asv,20)},striped=TRUE); output$taxa<-renderTable({req(rv$taxa);head(rv$taxa,20)},striped=TRUE)
  output$files<-renderUI({req(rv$run_dir);tags$ul(lapply(list.files(rv$run_dir,recursive=TRUE),tags$li))})
  output$log<-renderText(if(is.null(rv$log)||!file.exists(rv$log))"No run log yet." else paste(readLines(rv$log,warn=FALSE),collapse="\n"))

  observeEvent(input$run,{
    req(input$forward,input$training)
    if(input$mode=="paired") req(input$reverse)
    if(isTRUE(input$use_fastp)&&!exe_available("fastp")) return(showNotification("fastp was requested but is not installed or not on PATH.",type="error",duration=NULL))
    if(isTRUE(input$use_cutadapt)&&!exe_available("cutadapt")) return(showNotification("cutadapt was requested but is not installed or not on PATH.",type="error",duration=NULL))
    if(isTRUE(input$use_cutadapt)&&!nzchar(trimws(input$primer_f))) return(showNotification("Enter a forward primer for cutadapt.",type="error"))

    started<-Sys.time(); safe<-gsub("[^A-Za-z0-9_-]+","_",trimws(input$run_name)); if(!nzchar(safe))safe<-"FASTQ_Run"
    run_id<-paste0(format(started,"%Y%m%d_%H%M%S"),"_",safe); run_dir<-file.path(fastq_runs_root,run_id)
    for(d in c("01_QC","02_Preprocessed","03_Filtered","04_DADA2","05_ASV","06_Taxonomy","07_Export","Logs")) dir.create(file.path(run_dir,d),recursive=TRUE,showWarnings=FALSE)
    rv$run_dir<-normalizePath(run_dir,winslash="/",mustWork=TRUE); rv$log<-file.path(run_dir,"Logs","FASTQ_pipeline.log"); rv$status<-"Processing"
    err_obj<-NULL
    withProgress(message="Triple_A FASTQ",value=0,{tryCatch({
      fwd<-input$forward$datapath; fwd_names<-input$forward$name
      if(input$mode=="paired") { pairing<-pair_files(fwd_names,input$reverse$name); rev<-input$reverse$datapath[pairing$reverse_order]; samples<-pairing$samples } else { rev<-NULL; samples<-make.unique(sample_key(fwd_names)) }
      names(fwd)<-samples; if(!is.null(rev))names(rev)<-samples
      pre_f<-fwd; pre_r<-rev
      if(isTRUE(input$use_fastp)){
        incProgress(.08,detail="fastp"); out_f<-file.path(run_dir,"02_Preprocessed",paste0(samples,"_R1.fastq.gz")); out_r<-if(input$mode=="paired")file.path(run_dir,"02_Preprocessed",paste0(samples,"_R2.fastq.gz")) else NULL
        for(i in seq_along(samples)){ args<-c("-i",fwd[i],"-o",out_f[i],"--json",file.path(run_dir,"01_QC",paste0(samples[i],"_fastp.json")),"--html",file.path(run_dir,"01_QC",paste0(samples[i],"_fastp.html")),"--thread","2"); if(input$mode=="paired")args<-c(args,"-I",rev[i],"-O",out_r[i]); run_command("fastp",args,rv$log) }
        pre_f<-out_f; pre_r<-out_r
      }
      if(isTRUE(input$use_cutadapt)){
        incProgress(.08,detail="cutadapt"); out_f<-file.path(run_dir,"02_Preprocessed",paste0(samples,"_trim_R1.fastq.gz")); out_r<-if(input$mode=="paired")file.path(run_dir,"02_Preprocessed",paste0(samples,"_trim_R2.fastq.gz")) else NULL
        for(i in seq_along(samples)){ args<-c("-g",trimws(input$primer_f),"-o",out_f[i]); if(input$mode=="paired"){ if(nzchar(trimws(input$primer_r)))args<-c(args,"-G",trimws(input$primer_r)); args<-c(args,"-p",out_r[i],pre_f[i],pre_r[i]) } else args<-c(args,pre_f[i]); run_command("cutadapt",args,rv$log) }
        pre_f<-out_f; pre_r<-out_r
      }
      incProgress(.12,detail="DADA2 filtering")
      filt_f<-file.path(run_dir,"03_Filtered",paste0(samples,"_F.fastq.gz")); filt_r<-if(input$mode=="paired")file.path(run_dir,"03_Filtered",paste0(samples,"_R.fastq.gz")) else NULL
      if(input$mode=="paired") filtering<-dada2::filterAndTrim(pre_f,filt_f,pre_r,filt_r,truncLen=c(input$trunc_f,input$trunc_r),maxEE=c(input$max_ee_f,input$max_ee_r),minLen=input$min_len,truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=TRUE)
      else filtering<-dada2::filterAndTrim(pre_f,filt_f,truncLen=input$trunc_f,maxEE=input$max_ee_f,minLen=input$min_len,truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=TRUE)
      filtering<-as.data.frame(filtering); filtering$sample<-samples; write.csv(filtering,file.path(run_dir,"01_QC","filtering_summary.csv"),row.names=FALSE)
      incProgress(.18,detail="Learning errors"); errF<-dada2::learnErrors(filt_f,multithread=TRUE); saveRDS(errF,file.path(run_dir,"04_DADA2","error_forward.rds"))
      if(input$mode=="paired"){errR<-dada2::learnErrors(filt_r,multithread=TRUE);saveRDS(errR,file.path(run_dir,"04_DADA2","error_reverse.rds"))}
      incProgress(.22,detail="Inferring ASVs"); derepF<-dada2::derepFastq(filt_f);names(derepF)<-samples; ddF<-dada2::dada(derepF,err=errF,pool=if(input$pool)"pseudo" else FALSE,multithread=TRUE)
      if(input$mode=="paired"){derepR<-dada2::derepFastq(filt_r);names(derepR)<-samples;ddR<-dada2::dada(derepR,err=errR,pool=if(input$pool)"pseudo" else FALSE,multithread=TRUE);mergers<-dada2::mergePairs(ddF,derepF,ddR,derepR,minOverlap=input$min_overlap);seqraw<-dada2::makeSequenceTable(mergers);denoised<-vapply(mergers,function(x)sum(x$abundance),numeric(1))} else {seqraw<-dada2::makeSequenceTable(ddF);denoised<-vapply(ddF,function(x)sum(dada2::getUniques(x)),numeric(1))}
      seqtab<-dada2::removeBimeraDenovo(seqraw,method="consensus",multithread=TRUE);saveRDS(seqtab,file.path(run_dir,"04_DADA2","sequence_table_non_chimeric.rds"))
      incProgress(.18,detail="Assigning taxonomy"); taxa<-dada2::assignTaxonomy(seqtab,input$training$datapath,multithread=TRUE);saveRDS(taxa,file.path(run_dir,"06_Taxonomy","taxonomy_assignments.rds"))
      ids<-paste0("ASV",seq_len(ncol(seqtab))); abundance<-data.frame(ASV_ID=ids,t(seqtab),check.names=FALSE); write.csv(abundance,file.path(run_dir,"05_ASV","ASV_abundance_table.csv"),row.names=FALSE); write_fasta(colnames(seqtab),file.path(run_dir,"05_ASV","representative_sequences.fasta"))
      taxdf<-data.frame(ASV_ID=ids,as.data.frame(taxa,check.names=FALSE),check.names=FALSE);write.csv(taxdf,file.path(run_dir,"06_Taxonomy","ASV_taxonomy.csv"),row.names=FALSE,na="")
      prefixes<-c("k__","p__","c__","o__","f__","g__","s__"); lineage<-apply(taxa,1,function(x)paste(paste0(prefixes[seq_along(x)],x)[!is.na(x)&nzchar(x)],collapse=";")); triple<-data.frame(ASV_ID=ids,Taxonomy=lineage,t(seqtab),check.names=FALSE);write.csv(triple,file.path(run_dir,"07_Export","Triple_A_input_table.csv"),row.names=FALSE,na="")
      nonchim<-rowSums(seqtab); track<-data.frame(sample=samples,input=filtering[,1],filtered=filtering[,2],denoised=as.numeric(denoised[samples]),non_chimeric=as.numeric(nonchim[samples]),check.names=FALSE);write.csv(track,file.path(run_dir,"01_QC","read_tracking.csv"),row.names=FALSE)
      elapsed<-as.numeric(difftime(Sys.time(),started,units="secs")); summary<-data.frame(metric=c("Layout","Samples","Input reads","Filtered reads","Non-chimeric reads","ASVs","Elapsed"),value=c(input$mode,nrow(track),sum(track$input),sum(track$filtered),sum(track$non_chimeric),ncol(seqtab),fmt_time(elapsed)))
      metadata<-list(schema_version="5.0",module="FASTQ",purpose="preprocessing_only",status="completed",layout=input$mode,started=format(started,tz="UTC",usetz=TRUE),elapsed_seconds=elapsed,input_files=c(input$forward$name,input$reverse$name %||% character()),external_tools=list(fastp=isTRUE(input$use_fastp),cutadapt=isTRUE(input$use_cutadapt)),parameters=list(trunc_f=input$trunc_f,trunc_r=if(input$mode=="paired")input$trunc_r else NULL,max_ee_f=input$max_ee_f,max_ee_r=if(input$mode=="paired")input$max_ee_r else NULL,min_overlap=input$min_overlap,pool=input$pool),triple_a_input="07_Export/Triple_A_input_table.csv")
      jsonlite::write_json(metadata,file.path(run_dir,"run_metadata.json"),pretty=TRUE,auto_unbox=TRUE,null="null");saveRDS(list(table=triple,taxonomy=taxdf,tracking=track,metadata=metadata),file.path(run_dir,"07_Export","Triple_A_handoff.rds"))
      rv$tracking<-track;rv$asv<-triple;rv$taxa<-taxdf;rv$summary<-summary;rv$status<-paste("Completed in",fmt_time(elapsed));log_line(rv$log,rv$status);incProgress(.12,detail="Complete")
    },error=function(e)err_obj<<-e)})
    if(!is.null(err_obj)){rv$status<-paste("Failed:",conditionMessage(err_obj));log_line(rv$log,paste("ERROR",conditionMessage(err_obj)));showNotification(conditionMessage(err_obj),type="error",duration=NULL)}
  })
}
shinyApp(ui,server)
