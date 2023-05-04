library(phyloseq); packageVersion("phyloseq")
library(tidyverse); 
library(pheatmap)


## read the combined outpt file
emu <- read.delim("data/emu-combined-tax_id-counts.tsv")

emu %>% select(starts_with("barcode")) %>%
    round() %>% replace(is.na(.), 0) ->
    emuCounts 

emu %>% select(-starts_with("barcode")) %>%
    mutate_all(as.character) %>%
    mutate_all(~gsub("^$", NA, . )) %>%
    rev() ->
    emuTax

## ## Checking taxonomic anotation
table(emuTax$phylum)[order(table(emuTax$phylum))]
## obviously many fungi buat also Apicomplexa and Nematodes nicely
## covered
table(emuTax[emuTax$phylum%in%"Nematoda", "genus"])
## interesting: many credibly
table(emuTax[emuTax$phylum%in%"Chordata", "species"])
## complete rubbish: Why?
table(emuTax[emuTax$phylum%in%"Apicomplexa", "genus"])
## interesting: as usual many Gregarines, but also credible Eimeria
## and two Babesia
table(emuTax[emuTax$phylum%in%"Platyhelminthes", "genus"])
## interesting
table(emuTax[emuTax$phylum%in%"Streptophyta", "genus"])

## Now hte sample data
animals <- read.csv("data/AnimalData28S.csv")

animals$Barcode <- gsub("BC(\\d\\d) *", "barcode\\1", animals$Barcode)

rownames(animals) <- animals$Barcode

table(colnames(emuCounts)%in%rownames(animals))
## 4 sequencing barcodes detected that are not in our animals table

colSums(emuCounts)[colnames(emuCounts)%in%rownames(animals)]
colSums(emuCounts)[!colnames(emuCounts)%in%rownames(animals)]
## those have very, very few reads: NICE!

table(rownames(animals)%in%colnames(emuCounts))
## all samples in the animals table got sequences

## we can keep the counts for the proper samples only!
emuCounts <- emuCounts[colnames(emuCounts)%in%rownames(animals)]

## and bring the animals ("sample data") in the same order
animals <- animals[colnames(emuCounts), ]

PS <- phyloseq(otu_table(emuCounts, taxa_are_rows=TRUE ),
               tax_table(as.matrix(emuTax)),
               sample_data(animals))

## select PROJECTS, this should become independent scripts loading the
## PS objects as we progress with analyses

table(sample_data(PS)[, "project"])


Dupl <- subset_samples(PS, project%in%"Düppel")
Dupl <- subset_taxa(Dupl, taxa_sums(Dupl) > 1)

Dupl
## 203 taxa

NemDupl <- subset_taxa(Dupl, phylum%in%"Nematoda")
## 13 Nematodes 

pheatmap(log10(otu_table(NemDupl)+1),
         labels_row=tax_table(NemDupl)[, "genus"],
         labels_col=sample_data(NemDupl)$sample_ID,
         display_numbers=TRUE)


ApiDupl <- subset_taxa(Dupl, phylum%in%"Apicomplexa")
ApiDupl
## 4 Apicomplexa 

pheatmap(log10(otu_table(ApiDupl)+1),
         labels_row=tax_table(ApiDupl)[, "genus"],
         labels_col=sample_data(ApiDupl)$sample_ID,
         display_numbers=TRUE)


subset_taxa(PS, phylum%in%"Nematoda")
## of a total of 35, looks like flotation worked to an extent in the
## Düppel samples.

## Let's look at Lilla's cheetah blood and ticks. 

Lilla <- subset_samples(PS, project%in%"Lilla")
ApiLilla  <- subset_taxa(Lilla, phylum%in%"Apicomplexa")

ApiLilla <- subset_taxa(ApiLilla, taxa_sums(ApiLilla) >1)

pheatmap(log10(otu_table(ApiLilla)+1),
         labels_row=tax_table(ApiLilla)[, "genus"],
         labels_col=sample_data(ApiLilla)$sample_ID,
         display_numbers=TRUE)

### shit there's nothing!?
### In which semples are the interesting Apciomplexans

ApiPS  <- subset_taxa(PS, phylum%in%"Apicomplexa")
ApiPS <- subset_samples(ApiPS, colSums(otu_table(ApiPS)) >1)

ApiPS
## almost all (50/56) samples have Apicomplexans?! Most pobably only
## Gregarines

apply(tax_table(ApiPS), 2, table)

## "interesting Apicomplexans" 
ApiPS <- subset_taxa(PS, order%in%c("Eucoccidiorida", "Piroplasmida"))
ApiPS <- subset_samples(ApiPS, colSums(otu_table(ApiPS)) >1)

ApiPS
## Only 17 samples with "interesting" Apicomplexans
sample_data(ApiPS)
### All are rodents (but one Düppel ungulate). Brilliant! As expected
### (when missing out on blood parasites for Lilla)!


### What about Grace's Zebra samples, any Strongyles or Nematodes
### generally?

Grace <- subset_samples(PS, project%in%"Grace")
GraceNem <- subset_taxa(Grace, phylum%in%"Nematoda")
GraceNem <- subset_taxa(GraceNem, taxa_sums(GraceNem) > 1)

pheatmap(log10(otu_table(GraceNem)+1),
         labels_row=tax_table(GraceNem)[, "genus"],
         labels_col=sample_data(GraceNem)$sample_ID,
         display_numbers=TRUE)

tax_table(GraceNem)[, "species"]
## super nice!!
### Chabaudstrongylus ninhae" and "Haemonchus sp. 92619" are database
### artifacts but all these sequences are likely from strongyles!!!


### Finally Susana Soares' colourful mix of species, do we see
### something resembling the patterns for Düppel (flotations) and
### Grace (non-flotated Zebras)?

Sus <- subset_samples(PS, project%in%"SuzanaSuares")
SusNem <- subset_taxa(Sus, phylum%in%"Nematoda")
SusNem <- subset_taxa(SusNem, taxa_sums(SusNem) > 1)


png("figures/Grace_1st_heat.png")
pheatmap(log10(otu_table(SusNem)+1),
         labels_row=tax_table(SusNem)[, "genus"],
         labels_col=sample_data(SusNem)$sample_ID,
         display_numbers=TRUE)
dev.off()

## Yes, Petrovinema and Cylicocyclus again in Susana's Zebra. Also
## Haemonchus and (new here) Uncinaria in this Zebra. 

## Just for completeness the Apicomplexans for Susana:

SusApi <- subset_taxa(Sus, phylum%in%"Apicomplexa")
SusApi <- subset_taxa(SusApi, taxa_sums(SusApi) > 1)

pheatmap(log10(otu_table(SusApi)+1),
         labels_row=tax_table(SusApi)[, "genus"],
         labels_col=sample_data(SusApi)$sample_ID,
         display_numbers=TRUE)
## Only uninteresting Gregarines!
