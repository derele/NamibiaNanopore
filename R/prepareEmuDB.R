### packages we use
library(Biostrings); packageVersion("Biostrings")
library(taxonomizr); packageVersion("taxonomizr")
library(data.table)


## Downloading a very comprehensive rRNA database for
## [emu]("https://gitlab.com/treangenlab/emu") read classification from
## [RNACentral]("https://rnacentral.org/")'s ([ftp download "Active sequences are present in at least one expert database"][http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz)


## the file inconsistently uses T and U for DNA or RNA Thymin/Uracil
## bases we convert it to U only
## zcat rnacentral_active.fasta.gz | sed '/^[^>]/s/t/u/gi' > rnacentral_active_TU.fasta

## and then read it using Biostrings as an RNA string set
rRNA <- readRNAStringSet("/SAN/db/emuDB/ALLRna/rnacentral_active_TU.fasta")

### only rRNA
rRNA <- rRNA[grep("rRNA", names(rRNA))]

## for now using sequences longer than 400nt (have to check efficiency
## scaling of emu/minimap with size of databases to see if we have to
## remove more)
rRNA <- rRNA[width(rRNA)>400]

seqIDs <- gsub(" .*", "", names(rRNA))

## we also downloaded information on the RnaCentral ids
IDs <- fread("/SAN/db/emuDB/ALLRna/id_mapping.tsv.gz", sep="\t", 
             col.names=c("ID", "db", "ranges", "taxID", "type", "X"))
             
IDs <- IDs[IDs$ID%in%seqIDs, ]

tTax <- IDs[, .(.N), by = .(ID, taxID)]
## some IDs have multiple taxIDs
ttTax <- tTax[, .(.N), by = .(ID)]

## table(ttTax$N >1)
## ##    FALSE     TRUE 
## ## 13844855   141148 

## ## I guess we can forget about those for now
seqIDsWuniqTaxID <- ttTax[N==1, ID]

## alternatively (with a bit of coding effort) we could annotate those
## at deeper nodes of the taxonomy (last common ancestor).

## how many sequences are we losing with this?

## ## table(seqIDs%in%seqIDsWuniqTaxID)
##    FALSE     TRUE 
##   149677 13844855 

## nice, only 149677 sequences lost, that's okay for now.

rRNA <- rRNA[seqIDs%in%seqIDsWuniqTaxID]

## check whether the sequences are unique

## table(duplicated(rRNA))
## ## yes, they are

### align seqIDs and rRNA names again
seqIDs <- gsub(" .*", "", names(rRNA))

## keep only the IDs we need
setkey(IDs, ID)
IDs <- IDs[J(seqIDs)]

## to avoid confusion
rm(seqIDs)

## remove some remaining duplicates 
IDs <- unique(IDs[, list(ID, taxID)])

### Now we create an emu database from this:

## species_taxid.fasta database
##  sequences where each sequence header starts with the respective
##  species-level tax id (or lowest level above species-level if
##  missing) preceeding a colon [<species_taxid>:<remainder of
##  header>]

IDs[ , new := do.call(paste, c(list(taxID, ID), sep = ":"))]

## ## bring rRNA in the same order
## rRNA <- rRNA[order(IDs$ID)]


## if(any(gsub(" .*", "", names(rRNA)) != IDs$ID)){
##     stop("wrong order")
## }

## ## This doesn't work emu build database is needed
## names(rRNA) <-  IDs$new

## and the input fasta headers for that are simply the sequence IDs
names(rRNA) <-  IDs$ID

## Eukaryote subsetting
rRNA <- rRNA[names(rRNA)%in%IDs$ID]

## write it as DNA (with T instead of U) into a file
writeXStringSet(DNAStringSet(rRNA), "/SAN/db/emuDB/ALLRna/species_taxid.fasta")

## taxonomy.tsv tab separated datasheet of database taxonomy lineages
##  containing at columns: 'tax_id' and any taxonomic ranks
##  (i.e. species, genus, etc)
TaxDB <- "/SAN/db/taxonomy/update20200827/accessionTaxa.sql"

uTaxID <- unique(IDs[, taxID])

fullTax <- getTaxonomy(uTaxID, TaxDB)

## table(rowSums(apply(fullTax, 2, is.na)))

## ## 537594 entries have no NA at any level

## table(rowSums(apply(fullTax, 2, is.na))>3)
## ##  FALSE   TRUE 
## ## 618224  95342

## table(rowSums(apply(fullTax, 2, is.na))>3)
## ## 95342 have more than 3 levels NA... those are

badSpec <- table(fullTax[rowSums(apply(fullTax, 3, is.na))>1, "species"])
tail(badSpec[order(badSpec)])
## those are nothing special 

badSpec <- table(fullTax[rowSums(apply(fullTax, 2, is.na))>1, "species"])
tail(badSpec[order(badSpec)])
## also those are nothing special

fullTax <- as.data.frame(fullTax)

fullTax <- fullTax[fullTax$superkingdom%in%"Eukaryota",]

fullTax$tax_id <- gsub(" *", "", rownames(fullTax))

## as the database gets otherwise too big (and minimap2 then doesn't
## write @SQ headers into .sam) only select Eukaryota for now

write.table(fullTax[, c("tax_id", "superkingdom", "phylum", "class",
                        "order", "family", "genus", "species")],
            row.names=FALSE, "/SAN/db/emuDB/ALLRna/taxonomy.tsv",
            quote=FALSE, sep="\t")

## ## Now we still need seq2tax.map.tsv:  headerless two column
## ## tab-separated values, where each row contains (1) sequence
## ## header in database.fasta and (2) corresponding tax id"

## again the Eukaryote subsetting
IDs <- IDs[taxID%in%fullTax$tax_id, ]

write.table(IDs[, list(ID, taxID)], "/SAN/db/emuDB/ALLRna/seq2tax.map.tsv",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
            
### we can now create the database with the shell command: emu
### build-database ALLRna --sequences species_taxid.fasta --seq2tax
### seq2tax.map.tsv --taxonomy-list taxonomy.tsv"

### a call to emu could now look like this (re the database):
## emu abundance --db /SAN/db/emuDB/ALLRna/ALLRna/ 


## the database is now still too large, this means that the sam-file
## is written without @SQ header. We can make emu ignore this by
## changing EVERY OCCURENCE  of

##    sam_pysam = pysam.AlignmentFile(path)
## to
##    sam_pysam = pysam.AlignmentFile(path, check_sq=False)

## in the emu python script

### running into another error then back to above and shrink the
### database, as SQ lines should appear in sam once the db fits in
### memory
