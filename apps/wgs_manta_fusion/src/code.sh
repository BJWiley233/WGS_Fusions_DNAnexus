#!/bin/bash
# cram2bam 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://documentation.dnanexus.com/developer for tutorials on how
# to modify this file.


# dx run applet-GPkzVQ8JP6z8K131jjFbJP0F -y -icandidateSV_vcf=file-G9898X8JYGxjZVZY3x12F20F -icandidateSV_idx=file-G9898Y0JYGxpJ2fF5Vz0QJF6 -isubset_filtered_bam=file-G989580JFX7pJjKy4j162b5G -isubset_filtered_bai=file-G9895B8JFX7YGg2p4jJ7Gvgp  --destination project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Outputs/manta/full_workflow/ -iResults=Result.sai.mainscript.output.txt -iBackground=Background.sai.mainscript.output.txt

main() {
    set -ex -o pipefail
    # extract_fusion_matches_from_manta_vcf.pl
    echo "Value of perl_exe: '$perl_exe'"
    # chromoseq_translocations
    echo "Value of translocations: '$translocations'"
    echo "Value of genes_to_remove: '$genes_to_remove'"
    echo "Value of candidateSV_vcf: '$candidateSV_vcf'"
    echo "Value of candidateSV_idx: '$candidateSV_idx'"
    echo "Value of subset_filtered_bam: '$subset_filtered_bam'"
    echo "Value of subset_filtered_bai: '$subset_filtered_bai'"
    # results
    echo "Value of Result: '$Results'"
    echo "Value of Background: '$Background'"
    
    dx-download-all-inputs --parallel
    mv $candidateSV_idx_path /home/dnanexus/in/candidateSV_vcf
    mv $subset_filtered_bai_path /home/dnanexus/in/subset_filtered_bam

	
    perl $perl_exe_path $translocations_path $genes_to_remove_path $candidateSV_vcf_path $subset_filtered_bam_path $Results $Background


    res=$(dx upload $Results --brief)
    background=$(dx upload $Background --brief)
   
    dx-jobutil-add-output results "$res" --class=file
    dx-jobutil-add-output background "$background" --class=file
}