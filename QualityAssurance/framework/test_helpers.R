qa_find_root <- function(start=getwd()) { current<-normalizePath(start,winslash="/",mustWork=TRUE); repeat {if(file.exists(file.path(current,"Run_Triple_A.R"))&&file.exists(file.path(current,"AAApp","Common","Engine","Triple_A.R"))) return(current); parent<-dirname(current); if(identical(parent,current)) stop("Triple_A project root not found."); current<-parent}}
qa_expect_true <- function(x,message="Expectation was not TRUE") {if(!isTRUE(x)) stop(message,call.=FALSE); TRUE}
qa_expect_false <- function(x,message="Expectation was not FALSE") {if(!identical(x, FALSE)) stop(message,call.=FALSE); TRUE}
qa_expect_equal <- function(x,y,message="Values differ") {if(!isTRUE(all.equal(x,y))) stop(message,call.=FALSE); TRUE}
qa_expect_files <- function(paths) {missing<-paths[!file.exists(paths)&!dir.exists(paths)]; if(length(missing)) stop("Missing: ",paste(missing,collapse=", "),call.=FALSE); TRUE}
qa_write_reports <- function(root,results) {d<-file.path(root,"QualityAssurance","reports");dir.create(d,recursive=TRUE,showWarnings=FALSE);utils::write.csv(results,file.path(d,"validation_report.csv"),row.names=FALSE);saveRDS(results,file.path(d,"validation_report.rds"));jsonlite::write_json(results,file.path(d,"validation_report.json"),pretty=TRUE,na="null"); counts<-table(results$status); lines<-c("Triple_A Validation Report",paste("Generated:",Sys.time()),paste(names(counts),as.integer(counts),sep=": "),"",apply(results,1,function(x)paste0("[",x[['status']],"] ",x[['id']]," | ",x[['level']]," | ",x[['severity']]," | ",x[['description']],if(nzchar(x[['details']]))paste0(" | ",x[['details']])else"")));writeLines(lines,file.path(d,"validation_report.txt"));writeLines(capture.output(sessionInfo()),file.path(d,"session_info.txt"));invisible(results)}


qa_read_app_source <- function(root) {
  files <- c(
    file.path(root, "AAApp", "Biological", "app.R"),
    sort(list.files(
      file.path(root, "AAApp", "Biological", "modules"),
      pattern = "\\.[Rr]$",
      full.names = TRUE
    ))
  )
  files <- files[file.exists(files)]
  paste(
    unlist(lapply(files, readLines, warn = FALSE, encoding = "UTF-8")),
    collapse = "\n"
  )
}
