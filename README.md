# Analysis of Nanopore amplicon sequencing data generated during the GRK2046 field trip in Namibia


## Bioinformatic analysis and taxonomic annotation

We first ran guppy barcoder to assign the sequences to sample-specific
barcodes.

``` guppy_barcoder -i allPassedForNow -s Barcoded --barcode_kits
EXP-PBC096 --num_barcoding_threads 8 --fastq_out -q 0 ```

This produces a fastq file for each barcode for each of three
sequencing runs in the folder "Barcoded".

We put barcoded sequences from the diffrent sequencing runs into one
file

```for dir in $(ls -d barcode*);do cat $dir/*.fastq > $dir.fastq;
done```

We then use [emu]("https://gitlab.com/treangenlab/emu") for read
classification. Classification of Eukaryote sequences needs a special
database which we downloaded from
[RNACentral]("https://rnacentral.org/") see R/prepareEmuDB.R. The
focus of this database is on comprehensiveness.



```parallel -j 4 emu abundance --db /SAN/db/emuDB/ALLRna/ALLRna/
--output-dir AbuVSALLEuk --output-unclassified --keep-counts
--keep-files --threads 2 ::: Barcoded/*.fastq```

Some of this is very memory hungry and we therefore don't run more
than 2 threads on 4 files in parallel.


We then have to remove files with very few reads or assignments to
species. Those are actually unused barcodes that were mis-assigned
some (few) sequences.

```find . -type f -exec awk -v x=10 'NR==x{exit 1}' {} \; -exec rm -f
{} \;```

Then we can run another emu command to combine the abundance tables of
different barcodes.


