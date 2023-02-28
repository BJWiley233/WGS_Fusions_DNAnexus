# WGS_Fusions_DNAnexus

Contains 3 steps.

1. Samtools step to subset the bam file based on a fusion bed file such as in the resources folder of the `samtools_sv` app
Example
```bash
$ head fusion_genes.bed
chr1    2228694 2310119
chr1    3069210 3438621
chr1    6785324 7769706
chr1    9234775 9271337
chr1    9292880 9369532
chr1    10398592        10420144
chr1    15756607        15786594
chr1    15847864        15940460
chr1    22001656        22012539
chr1    26693236        26782110
```

2. 



### Docker Images required are upload as tar files in the `dockers` folder and can be unloaded.  Just change the `$dnanexus_link` in the `dxapp.json` files for upload tar files to DNAnexus
