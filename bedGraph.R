works_with_R("3.1.2",
             dplyr="0.4.0",
             "Rdatatable/data.table@200b5b40dd3b05112688c3a9ca2dd41319c2bbae")

colClasses <- c("chrom"="factor",
                "bases"="integer",
                file="NULL")
chroms <-
  read.table("~/projects/seq-segment/chromInfo.txt.gz",
             colClasses=colClasses,
             nrow=24,
             row.names=1,
             col.names=names(colClasses))

colClasses <-
  c("snp"="NULL",
    "chr"="factor",
    "position"="integer",
    "logratio"="numeric",
    "mean"="NULL",
    "Bkp"="NULL",
    "Gnl"="NULL",
    "Out"="NULL",
    foo="NULL")
NB1388 <- fread("NB1388T_D.Copy_Number_snp.csv",
                colClasses=as.character(colClasses))
setnames(NB1388, names(colClasses)[colClasses != "NULL"])
NB1388$AD <- fread("V_AD_NB1388T_D.txt")
NB1388$LRR <- fread("V_LRR_NB1388T_D.txt")
with(NB1388, stopifnot(logratio == LRR))

header.tmp <-
  paste('track',
        'type=bedGraph',
        'db=hg19',
        'export=yes',
        'visibility=full',
        'maxSegments=20',
        'alwaysZero=on',
        'share=curie.fr',
        'graphType=points',
        'yLineMark=0',
        'yLineOnOff=on',
        'name=%s%s',
        'description="%s %s"')
nb.id <- "NB1388"

some.chr <- NB1388 %>%
  mutate(chrom=factor(paste0("chr", chr),
           paste0("chr", c(1:22, "X", "Y")))) %>%
  filter(!is.na(chrom)) %>%
  mutate(bases=chroms[as.character(chrom), "bases"],
         end.int=ifelse(position > bases, bases, position),
         start.int=end.int - 1L,
         end.chr=sprintf("%d", end.int),
         start.chr=sprintf("%d", start.int))
subset(some.chr, position > bases) # no positions bigger than hg19 chrom ends.

for(data.type in c("AD", "LRR")){
  header <- sprintf(header.tmp, nb.id, data.type, nb.id, data.type)
  file.name <- paste0(nb.id, data.type, ".bedGraph.gz")
  con <- gzfile(file.name, "w")
  writeLines(header, con)
  df <- data.frame(some.chr)[, c("chrom", "start.chr", "end.chr", data.type)]
  write.table(df, con, quote=FALSE, row.names=FALSE, col.names=FALSE)
  close(con)
}
