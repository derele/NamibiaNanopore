 Analysis of Nanopore amplicon sequencing data generated during the GRK2046 field trip in Namibia
====================================================================================

You might want to skip the bioinformatics section, which is mostly for
Marly, helping her to implement similar analyses for other MinIon data
but also might later - in less detail - appear in methods parts.

Project specific sections are: 
- [Grace's project](##graces-project)
- ["Düppel" project](##düppel-project)
- [Susana's project](##susanas-project)
- [Lilla's project](##lillas-project)

## Bioinformatic analysis and taxonomic annotation

We first ran guppy barcoder to assign the sequences to sample-specific
barcodes.

``` 
guppy_barcoder -i allPassedForNow -s Barcoded --barcode_kits
EXP-PBC096 --num_barcoding_threads 8 --fastq_out -q 0 
```

This produces a fastq file for each barcode for each of three
sequencing runs in the folder "Barcoded".

We put barcoded sequences from the diffrent sequencing runs into one
file (for each barcode):

```
for dir in $(ls -d barcode*);do cat $dir/*.fastq > $dir.fastq;
done
```

We then use [emu]("https://gitlab.com/treangenlab/emu") for read
classification. Classification of Eukaryote sequences needs a special
database which we downloaded from
[RNACentral]("https://rnacentral.org/") see
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
were identified and their respective coutns in the processed
samples. The taxonomic identification (at this stage) is very
dependent on the content of the database: worms like Haemonchus
contortus, which have the complete rDNA region
(18S-ITS1-5.8S-ITS2-28S) we sequenced in the database will "attract" a
lot of the annotations. This obviously doesn't mean that we detected
Haemonchus in a Zebra, but rather that there's something most similar
to it for the amplified region. 

In the next stage of the analysis I'll extract all the sequences for
the interesting taxa (by this first annotation) and produce alignemnts
and consensus sequences for them. Interested collaborators in their
individual projects can then "manually" revise the taxonomic
annotation for these sequences, or e.g. build phylogentic trees from
them (e.g. for the diverse small strongyles).


The following project sepecific plots ("heatmaps") show read numbers
(log10(x+1) transformed) for samples (IDs you used when you
"submitted" the samples for PCR and sequencing) in the colums and taxa
in the rows.  


## Grace's project
### detection of (small strongyle) nematodes in crude Zebra feces

![Log10 +1 transformed counts for
nematodes](https://github.com/derele/NamibiaNanopore/blob/main/figures/Grace_1st_heat.png)

"Chabaudstrongylus ninhae" and "Haemonchus sp. 92619" are database
artifacts (those are likely other Strongyles, but these species are
not covered well enough in the database) but all these sequences are
likely from strongyles!!! We will see in alignemtns how many credible
sequence variants these four nematode taxa correspond to. 


## Düppel-project (Zaneta, Luis, Joshua and Philipp)
### detection of nematodes in ungulate fecal flotations


![Log10 +1 transformed counts for
nematodes](https://github.com/derele/NamibiaNanopore/blob/main/figures/DuplNem_1st_heat.png)

In total out of 203 taxa detected in those samples (mostly fungi, some
plants). Detecting a total of 11392 for 13 Nematode taxa in this is an
enrichment compared to other fecal samples (see below).

If I remember correclty also some (few, two?) of the "Düppel samples"
(which IDs?)  were not flotated but crude fecal mater, correct? So we
could use this in a first analysis of how flotation enriches for the
PCR products of Nematodes.


## Susana's-project
### detection of parasites in curede fecal samples of multiple species

Susana Soares' colourful mix of species, allows us to ask a specific
questions related to other projects: 

![Log10 +1 transformed counts for
nematodes](https://github.com/derele/NamibiaNanopore/blob/main/figures/Susana_1st_heat.png)

Do we see the same in the (one) Zebra as in Grace's Zebra?

Yes, Petrovinema and Cylicocyclus again in Susana's Zebra (SS7Z is a
Zebra). Also Haemonchus and (new here) Uncinaria in this Zebra.

Do we see something resembling the patterns of flotated materials from
ungulates (the Düppel project) in ungulate flotations of crude feces?
We see much fewer Nematodes in these samples. Maybe the lack of
flotation or the different host species? Cooperia and Haemonchus in
the Ostrich (SS3OS) and Haemonchus in the Giraffe (SS6G) and in SS2Ro
are interesting!?


## Lilla's-project
### detection of parasites in ticks and cheetah blood
