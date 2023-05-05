 Analysis of Nanopore amplicon sequencing data generated during the GRK2046 field trip in Namibia
====================================================================================

You might want to skip the bioinformatics section, which is mostly for
Marly, helping her to implement similar analyses for other MinIon data
but also might later - in less detail - appear in methods parts.
Project specific sections are: 
- [Grace's project](#graces-project)
- ["Düppel" project](#düppel-project)
- [Susana's project](#susanas-project)
- [Lilla's project](#lillas-project)
- [Mouse project](#mouse-project)

## Bioinformatic analysis and taxonomic annotation

We first ran guppy barcoder to assign the sequences to sample-specific
barcodes.

``` 
guppy_barcoder -i allPassedForNow -s Barcoded --barcode_kits
EXP-PBC096 --num_barcoding_threads 8 --fastq_out -q 0 
```

This produces a fastq file for each barcode for each of three
sequencing runs in the folder "Barcoded".

We put barcoded sequences from the different sequencing runs into one
file (for each barcode):

```
for dir in $(ls -d barcode*);do cat $dir/*.fastq > $dir.fastq;
done
```

We then use [emu](https://gitlab.com/treangenlab/emu) for read
classification. Classification of Eukaryote sequences needs a special
database which we downloaded from
[RNACentral](https://rnacentral.org/) see
[R/prepareEmuDB.R](https://github.com/derele/NamibiaNanopore/blob/main/R/prepareEmuDB.R). The
focus of this database is on comprehensiveness, in this case for all
rRNA of Eukaryotes.


```
parallel -j 4 emu abundance --db /SAN/db/emuDB/ALLRna/ALLRna/
--output-dir AbuVSALLEuk --output-unclassified --keep-counts
--keep-files --threads 2 ::: Barcoded/*.fastq
```

The minimap step of this is memory hungry and we therefore don't run
more than 2 threads on 4 files (8 threads total) in parallel.

We then have to remove files with very few reads or assignments to
species. Those are actually unused barcodes that were mis-assigned
some (few) sequences.

```
find . -type f -exec awk -v x=10 'NR==x{exit 1}' {} \; -exec rm -f
{} \;
```

Then we can run another emu command to combine the abundance tables of
different barcodes.

```
 emu combine-outputs --counts AbuVSALLEuk/ tax_id
 ```

This data is now in the '/data' folder here in this repository
(emu-combined-tax_id-counts.tsv)[https://github.com/derele/NamibiaNanopore/blob/main/data/emu-combined-tax_id-counts.tsv]. We
read it in the script
(emuResults.R)[https://github.com/derele/NamibiaNanopore/blob/main/R/emuResults.R],
to produce the analyses below.

# Individual projects

The following graphs give a first impression which database sequences
were identified and their respective counts in the processed
samples. The taxonomic identification (at this stage) is very
dependent on the content of the database: worms like *Haemonchus
contortus*, which have the complete rDNA region
(18S-ITS1-5.8S-ITS2-28S) we sequenced in the database will "attract" a
lot of the annotations. This obviously doesn't mean that we detected
*Haemonchus* in a Zebra, but rather that there's something most
similar to it for the amplified region.

In the next stage of the analysis I'll extract all the sequences for
the interesting taxa (by this first annotation) and produce alignments
and consensus sequences for them. Interested collaborators in their
individual projects can then "manually" revise the taxonomic
annotation for these sequences, or e.g. build phylogenetic trees from
them (e.g. for the diverse small strongyles).


The following project specific plots ("heatmaps") show read numbers
(log10(x+1) transformed) for samples (IDs you used when you
"submitted" the samples for PCR and sequencing) in the columns and taxa
in the rows.  


## Grace's project
### detection of (small strongyle) nematodes in crude Zebra feces

![Log10 +1 transformed counts for
nematodes](https://github.com/derele/NamibiaNanopore/blob/main/figures/Grace_1st_heat.png)

*Chabaudstrongylus ninhae* and *Haemonchus sp. 92619* are database
artifacts (those are likely other strongyles, but these species are
not covered well enough in the database) but all these sequences are
likely from strongyles!!! We will see in alignments how many credible
sequence variants these four nematode taxa correspond to.


## Düppel-project
###  (Zaneta, Luis, Joshua and Philipp) detection of nematodes in ungulate fecal flotations


![Log10 +1 transformed counts for
nematodes](https://github.com/derele/NamibiaNanopore/blob/main/figures/DuplNem_1st_heat.png)

In total out of 203 taxa detected in those samples (mostly fungi, some
plants). Detecting a total of 11392 for 13 Nematode taxa in this is an
enrichment compared to other fecal samples (see below).

If I remember correctly also some (few, two?) of the "Düppel samples"
(which IDs?)  were not flotated but crude fecal mater, correct? So we
could use this in a first analysis of how flotation enriches for the
PCR products of Nematodes.


## Susana's project
### detection of parasites in crude fecal samples of multiple species

Susana Soares' colorful mix of species, allows us to ask a specific
questions related to other projects:

![Log10 +1 transformed counts for
nematodes](https://github.com/derele/NamibiaNanopore/blob/main/figures/Susana_1st_heat.png)

Do we see the same in the (one) Zebra as in Grace's Zebra?

Yes, *Petrovinema* and *Cylicocyclus* again in Susana's Zebra (SS7Z is
a Zebra). Also *Haemonchus* and (new here) *Uncinaria* in this Zebra.

Do we see something resembling the patterns of flotated materials from
ungulates (the Düppel project) in ungulate flotations of crude feces?
We see much fewer Nematodes in these samples. Maybe the lack of
flotation or the different host species? Cooperia and *Haemonchus* in
the Ostrich (SS3OS) and *Haemonchus* in the Giraffe (SS6G) and in SS2Ro
are interesting!?


## Lilla's project
### detection of parasites in ticks and cheetah fecal samples

So far only seqeunces annotated to the genus "Gregarina" are found
those are database artefacts (or, less likely imo, PCR artefacts). We
have to investigate some of these sequences more thoroughly, but for
now: we should establish more specific Apicomplexa primers to really
amplify from blood or ticks.

Regarding the latter, we found many tick related arthropod sequences.

![Log10 +1 transformed counts for
arthropods](https://github.com/derele/NamibiaNanopore/blob/main/figures/Lilla_Arthro_1st_heat.png)

The species in the hyena feces are credible diet or sampling
contamination. The other arthropods are all tick sequences. As this
should be from a single species it will be helpful to establish
alignment clustering procedures to derive a consensus sequence (in
this case one).

## Mouse project
### Fay's,  Marly's and Otto's mouse project

![Log10 +1 transformed counts for
parasites](https://github.com/derele/NamibiaNanopore/blob/main/figures/Rodents_1st_heat.png)

This figure nicely shows the challenges that still persist: many
detections seem to be split over different taxonomic annotations. This
is likely asa we see closely related taxa appear in the same sample.
The true species was not available in the database and the emu
algorithm assigned them to different (similarly distant) taxa.
Therefor we need to extract the sequences (spererately for matches to
reference at certain taxonomic levels), and then align and re-cluster
them without reference.

For the mouse project this also shows nicely that we likely have a
couple of different *Eimeria* species: four or five would be my guess
based on the abundance patterns alone.

