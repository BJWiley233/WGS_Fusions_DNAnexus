# WGS_Fusions_DNAnexus

Can detect fusions from low depth WGS data.  Written for DNAnexus as a workflow.

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

3. Update the perl execution script from WashU in the `wgs_manta_fusion` app `dxapp.json` file.  The script is located in apps/wgs_manta_fusion/resources.
Make sure to upload the `translocations` and `genes_to_remove` file also located apps/wgs_manta_fusion/resources.



#### Docker Images required are upload as tar files in the `dockers` folder and can be unloaded.  Just change the `$dnanexus_link` in the `dxapp.json` files for upload tar files to DNAnexus
