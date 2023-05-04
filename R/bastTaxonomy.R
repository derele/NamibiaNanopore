library(taxonomizr)

TaxDB <- "/SAN/db/taxonomy/update20200827/accessionTaxa.sql"

## BLAf <- "/SAN/MinKnowNamibia/ParasiteLongAmprDNA/BarCodedSamples/barcode13_300selectionVSnt.blt"

### a RODENT
BLAf <- "/SAN/MinKnowNamibia/ParasiteLongAmprDNA/BarCodedSamples/barcode04_1000sectionVSnt.blt"

### a Zebra
## BLAf <-  "/SAN/MinKnowNamibia/ParasiteLongAmprDNA/BarCodedSamples/barcode89_300sectionVSnt.blt"

BLA <- read.delim(BLAf, header=FALSE)

names(BLA) <-  c("qaccver", "saccver", "pident", "length",
                 "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "staxid")

summary(BLA$length)

fullTax <- getTaxonomy(BLA$staxid, TaxDB)

BLATax <- cbind(BLA, fullTax)

BLATaxU <- BLATax[!duplicated(BLATax$qaccver), ]

table(BLATax[!duplicated(BLATax$qaccver), "phylum"])

table(BLATax[,"phylum"])

hist(BLATax[, "pident"])


tapply(BLATaxU$pident, BLATaxU$phylum, mean)

tapply(BLATaxU$length, BLATaxU$phylum, mean)
tapply(BLATaxU$length, BLATaxU$phylum, max)

