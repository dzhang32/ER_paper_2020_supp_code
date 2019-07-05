## Function to re-format the data from GTEx as well the validation from GTEx

load("/home/sguelfi/projects/R/hipp/data/expression/splitReads/splitReads.filtered.annotated.rda")

GTEx.info <- read.delim("/data/recount/GTEx/sample_ids.tsv",header = F)
recount.info <- read.delim("/data/recount/GTEx/SRP012682.tsv",header = T)

library(dplyr)

GTEx.info <- GTEx.info %>% 
  filter(V2=="SRP012682") %>%  
  mutate(V3_chr = as.character(V3)) %>%  
  left_join(recount.info %>% 
              mutate(run_chr = as.character(run)) %>% 
              select(run_chr,smts,smtsd),
            by=c("V3_chr"="run_chr"))


tissues <- names(table(GTEx.info$smtsd))

##### splits the junction coverage file by the sample IDs per tissue ready to be reformatted in later step
## select only those in parentesis
#tissues <- tissues[grep("(",x = tissues,fixed = T)]
sapply(tissues,function(x){
  
  write.table(GTEx.info %>%  filter(smtsd==x) %>% select(V1),file="/data/recount/GTEx/toGrep.txt",col.names =F,row.names = F)
  system(paste0("cat /data/recount/GTEx/SRP012682.junction_coverage.tsv | grep -w -f /data/recount/GTEx/toGrep.txt > /data/recount/GTEx/gtex_junction_coverage_by_tissue/",
                tissues <- gsub('\\)','',gsub('\\(','_',gsub(" ","",x))))) 
  message(paste(Sys.time(),":",x,"Done"))            
  
})


library(data.table)
library(stringr)
library(doMC)
library(itertools)

## get the coverage
# tissues <- gsub('\\(','_',names(table(GTEx.info$smtsd)))
# tissues <- gsub('\\)','',tissues)
## run this for those that have parentesis
#tissues <- tissues[grep("_",x = tissues,fixed = T)]
#idea taken from http://tomlogan.co.nz/blogs/benchmark_parallel_bigdataframe.html
sapply(tissues,function(x){
  
  message(paste(Sys.time(),x))
  message(paste(Sys.time(),"Loading table"))
  GTEx.info.tmp <- GTEx.info %>%  filter(smtsd==x) %>% select(V1)
  ## remove the parantesis
  x <- gsub('\\(','_',x)
  x <- gsub('\\)','',x)
  dir.create(paste0("/data/recount/GTEx/gtex_junction_coverage_by_tissue/datatable/",gsub(" ","",x)),showWarnings = F)
  tmp <- fread(paste0("/data/recount/GTEx/gtex_junction_coverage_by_tissue/",gsub(" ","",x)))
  message(paste(Sys.time(),"Updating table"))
  
  registerDoMC(15)
  pb <- txtProgressBar(min = 0, max = 15, style = 3)
  invisible({foreach(m=isplitRows(tmp, chunks=15),j=icount(),.combine='c') %dopar% {
    dt <- data.table(matrix(nrow=nrow(m),ncol=nrow(GTEx.info.tmp)+1))
    dt <- dt[, lapply(.SD, as.numeric)]
    setnames(dt,  c("juncID",as.character(GTEx.info.tmp$V1) ))
    for(i in 1:nrow(m))
    {
      i.tmp <- data.table(cbind(V2=as.character(unlist(str_split(m[i,"V2"], ","))),V3=as.numeric(unlist(str_split(m[i,"V3"], ","))))) %>% 
        filter(V2 %in% colnames(dt))
      dt[i, c("juncID",as.character(unlist(str_split(i.tmp[,"V2"], ",")))):= as.list(c(as.numeric(unlist(str_split(m[i,"V1"], ","))),as.numeric(str_split(i.tmp[,"V3"], ","))))]
      rm(i.tmp)
    }
    fwrite(dt,file=paste0("/data/recount/GTEx/byTissue/datatable/",gsub(" ","",x),"/tmp.",j))
    rm(dt)
    setTxtProgressBar(pb, j)
  }})
  message(paste(Sys.time(),"Merging files"))
  system(paste0("cat /data/recount/GTEx/byTissue/datatable/",gsub(" ","",x),"/tmp.1 | head -n1 > /data/recount/GTEx/byTissue/datatable/",gsub(" ","",x),"/",gsub(" ","",x),".csv"))
  system(paste0("for f in /data/recount/GTEx/byTissue/datatable/",gsub(" ","",x),"/tmp.*; do cat $f | tail -n +2 >> /data/recount/GTEx/byTissue/datatable/",gsub(" ","",x),"/",gsub(" ","",x),".csv; done"))
  rm(tmp,GTEx.info.tmp)
  message(paste(Sys.time(),"Task completed"))
})

## to merge the files  
# cat a-randomly-selected-csv-file.csv | head -n1 > merged.csv
# for f in *.csv; do cat "`pwd`/$f" | tail -n +2 >> merged.csv; done 

library(genomation)
## load junctions from GTEx, the link to the http link file can be use directly
GTEx.bed <- readBed("/data/recount/GTEx/SRP012682.junction_id_with_transcripts.bed.gz",zero.based = F)
load("~/projects/R/hipp/data/expression/splitReads/splitReads.filtered.annotated.rda")

## we add -1 to the coordinates to match the UCSC coordinates.
library(GenomicRanges)

tmp.table.GR <- GRanges(paste0("chr",as.character(tmp.table$chr)),
                        IRanges(start=tmp.table$start-1,tmp.table$stop-1),
                        strand = tmp.table$strand,
                        juncId=tmp.table$junId)

# tmpOverlap <- countOverlaps(tmp.table.GR, GTEx.bed,
#                             maxgap=0L, minoverlap=1L,
#                             type="equal",
#                             ignore.strand=FALSE)

tbOverlap <- findOverlaps(tmp.table.GR, GTEx.bed,
                          maxgap=0L, minoverlap=1L,
                          type="equal",
                          ignore.strand=FALSE)


library(stringr)

tbOverlap <- as.data.frame(tbOverlap)

head(tbOverlap)

idx.overlap <- cbind(tmp.table.GR[tbOverlap$queryHits]$juncId,
                     unlist(lapply(str_split(GTEx.bed[tbOverlap$subjectHits]$name,"\\|"),function(x){return(x[1])})))

library(data.table)

## load the validation table
tissues <- names(table(GTEx.info$smtsd))
tissues <- gsub('\\(','_',tissues)
tissues <- gsub('\\)','',tissues)
for (i in tissues)
{
  message(paste(Sys.time(),i))
  val.tab <- fread(paste0("/data/recount/GTEx/byTissue/datatable/",gsub(" ","",i),"/",gsub(" ","",i),".csv"),header = T)
  #length(intersect(idx.overlap,val.tab$juncID))/nrow(tmp.table)
  val.tab <- val.tab[apply(val.tab,1,function(x){
    sum(x[2:length(x)]>0,na.rm = T)> (length(x)-1)*0.05
  }),]
  tmp.table[which(as.character(tmp.table$junId) %in% (idx.overlap[which(idx.overlap[,2] %in%  as.character(val.tab$juncID)),1])) ,gsub(" ","",i)] <- 1  
  rm(val.tab)
}

head(GTEx.info)
head(tmp.table)


#save(tmp.table,file="~/projects/R/hipp/data/expression/splitReads/jun.table.ann.GTEx.rda")
load(file="~/projects/R/hipp/data/expression/splitReads/jun.table.ann.GTEx.rda")





