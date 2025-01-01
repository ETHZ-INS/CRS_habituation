library(epiwraps)
lf <- list.files("aligned", pattern="bam$", full=TRUE)
names(lf) <- gsub("\\.bam$","",basename(lf))
for(f in names(lf)){
  print(f);
  f2 <- paste0("tracks_ext3/",f,".bw")
  if(file.exists(f2)){
    message("Skipping ", f2)
  }else{
    bam2bw(lf[[f]], f2, paired=TRUE, binWidth=10L, shift=c(4L,-5L), type="ends", extend=3L, includeDuplicates=FALSE)
  }
}
