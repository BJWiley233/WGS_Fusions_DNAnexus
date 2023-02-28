# WGS_Fusions_DNAnexus

Can detect fusions from low depth WGS data using manta and scripts from Washington University in St. Louis.  Written for DNAnexus as a workflow.

<br>

Contains 3 steps.

1. Samtools step to subset the bam file based on a fusion bed file such as in the resources folder of the `samtools_sv` app

Example
```bash
$ head apps/samtools_sv/resources/fusion_genes.bed
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


2. Call manta structural variants with subset bam and manta config file located in apps/manta_svresources/configManta_1.py.ini.  Must upload config and edit in dxapp.json else you can use after resources is untarred to the DNAnexus instance.


3. Update the perl execution script from WashU in the `wgs_manta_fusion` app `dxapp.json` file.  The script is located in apps/wgs_manta_fusion/resources.
Make sure to upload the `translocations` and `genes_to_remove` file also located apps/wgs_manta_fusion/resources. Must upload and edit in dxapp.json else you can use after resources is untarred to the DNAnexus instance.

Example of translocations
```bash
$ head apps/wgs_manta_fusion/resources/chromoseq_translocations.withoutselfgenes_plus_solid_tumor.txt
chr1	3069210	3438621	chr1	2228694	2310119	PRDM16_SKI	.	+	+
chr1	110338505	110346677	chr22	40410280	40636685	RBM15_MRTFA	.	+	-
chr1	148808465	149032955	chr5	150113836	150155860	PDE4DIP_PDGFRB	.	+	-
chr1	154157317	154192058	chr5	150113836	150155860	TPM3_PDGFRB	.	-	-
chr1	186311651	186375325	chr8	38411138	38468834	TPR_FGFR1	.	-	-
chr1	221701423	221742176	chr1	3069210	3438621	DUSP10_PRDM16	.	-	+
chr1	234604268	234609525	chr17	40309193	40356796	IRF2BP2_RARA	.	-	+
chr2	32946971	33399359	chr2	32357027	32618899	LTBP1_BIRC6	.	+	+
chr2	43230835	43596046	chr3	169084760	169663470	THADA_MECOM	.	-	-
chr2	54456316	54671445	chr5	150113836	150155860	SPTBN1_PDGFRB	.	+	-
```


#### Docker Images required are upload as tar files in the `dockers` folder and can be unloaded.  Just change the `$dnanexus_link` in the `dxapp.json` files for upload tar files to DNAnexus

#### Most steps require a reference, reference index, and reference dictionary to be upload to DNAnexus and added to the inputs of the applet.