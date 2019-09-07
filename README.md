# Jose_2019

**Prerequisites:**  
[cutadapt 2.5](https://cutadapt.readthedocs.io/en/stable/index.html)  
[STAR-2.7.2b](https://github.com/alexdobin/STAR)  
[gffread utility](http://ccb.jhu.edu/software/stringtie/gff.shtml)  

Transcriptome (polyA captured mRNA-seq) samples were sequenced in paired-end 150 nt mode on Illumina sequencer.
Raw sequencing files are available from [GEO]().

### Preparing genome annotation and index files
Mouse genomic sequences and annotation files (GRCm38.p6) were downloaded from the [NCBI repository](ftp://ftp.ncbi.nih.gov/genomes/M_musculus/). To obtain genome assembly, download fasta files of individual chromosomes from ```Assembled_chromosomes/seq/``` folder and concatenate them in the ascending order (omit mitochondrial chromosome and sex chromosomes). Edit chromosome names to match annotation names in gff3 (for example convert ```>ref|NC_000067.6|``` to ```>NC_000067.6```) and say thanks to the NIH staff for letting you do it.  

| files               | MD5 check sum (unzipped)         | Description                                               |
| ------------------- |:--------------------------------:| ----------------------------------------------------------|
| GRCm38.p6.rna.fa    | b22037bd202465ce84cee9f6409331e8 | RNA in fasta format, coding + noncoding                   |
| GRCm38.p6.genome.fa | 49e0e5a638620d90990be88c81030923 | Genome sequence (nuclear genome only, no sex chromosomes) |
| GRCm38.p6.gbk       | adc1125bf6b2c3b5a52414fe2fe98ac7 | RNA in gene bank format, coding + noncoding               |
| GRCm38.p6.gff3      | ab982471b0b29ebde3d966ec2424253f | Genome annotation                                         | 


<details><summary><b>Customizing genome annotation</b></summary>  

**Customize genome annotation**  
Annotation of extrachromosomal contigs and sex chromosomes was omitted. 'Gnomon' (Predicted) records from gff file were also omitted and only 'RefSeq' and 'BestRefSeq' (manually curated) kept. Perl and R scripts are included in the GitHub repository.   
```bash
Discard_extrachromosomal_annotation.pl GRCm38.p6.gff3 >GRCm38.p6.custom.gff
Discard_gnomon_annotation.pl >GRCm38.p6.Refseq.gff	# automatically takes GRCm38.p6.custom.gff as an input
```
**Remove non-coding RNA genes**, leave only coding genes with their mRNA, transcript, exon, and CDS children. Fix the gff annotation from previous script by matching gene coordinates with the childern coordinates (occured due to removal of Gnomon features).
```bash
Discard_noncoding_annotation.R
```

**Convert annotation from GFF3 to GTF format**  
```bash
gffread GRCm38.p6.Refseq.coding.gff -T -o GRCm38.p6.Refseq.coding.gtf
# -T          - convert gff/gtf
```
</details>


<details><summary><b>Building genomic index</b></summary>
  
```bash  
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ./Mouse_index/ --genomeFastaFiles ./GRCm38.p6.genome.fa --sjdbGTFfile ./GRCm38.p6.Refseq.coding.gtf
```
</details>

