#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail

    echo ${bam_path}
    echo ${bam_index_path}
    echo ${manta_config}
    echo ${reference_path}
    echo ${reference_dict_path}
    echo ${reference_fai_path}
    echo ${dockerimage_manta_sv_path}
    echo ${non_wgs}
    echo ${manta_output_contigs}

    ls /usr/local/bin/manta_sv.sh

    dx-download-all-inputs --parallel

    mv $bam_index_path ~/in/bam
    mv $reference_path .
    mv $reference_fai_path .
    mv $reference_dict_path .
    mv $manta_config_path .

    echo $HOME
    sed -i "s|/illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa|${HOME}/${reference_name}|" $manta_config_name
    cat $manta_config_name



    docker load -i ${dockerimage_manta_sv_path}
    ## get <eid>_23153_0_0 from <eid>_23153_0_0.bqsr.bam
    eid_nameroot=$(echo $bam_name | cut -d'.' -f1)
    
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /usr/local/reference/:/usr/local/reference/ -w /home/dnanexus mgibio/manta_somatic-cwl:1.6.0 \
        /bin/bash /usr/local/bin/manta_sv.sh ${bam_path} ${HOME}/${reference_name} ${HOME}/${manta_config_name} ${non_wgs} ${manta_output_contigs}

    echo "results: "
    ls results/variants
    
    #File? diploid_variants = "results/variants/diploidSV.vcf.gz"
    # diploid_variants=$(dx upload "results/variants/diploidSV.vcf.gz" --brief)
    # dx-jobutil-add-output diploid_variants --class=file "$diploid_variants"

    #File? diploid_variants_tbi = "results/variants/diploidSV.vcf.gz.tbi"
    # diploid_variants_tbi=$(dx upload "results/variants/diploidSV.vcf.gz.tbi" --brief)
    # dx-jobutil-add-output diploid_variants_tbi --class=file "$diploid_variants_tbi"
    
    #File? somatic_variants = "results/variants/somaticSV.vcf.gz"
    # somatic_variants=$(dx upload "results/variants/somaticSV.vcf.gz" --brief)
    # dx-jobutil-add-output somatic_variants --class=file "$somatic_variants"

    #File? somatic_variants_tbi = "results/variants/somaticSV.vcf.gz.tbi"
    # somatic_variants_tbi=$(dx upload "results/variants/somaticSV.vcf.gz.tbi" --brief)
    # dx-jobutil-add-output somatic_variants_tbi --class=file "$somatic_variants_tbi"
    
    #File all_candidates = "results/variants/candidateSV.vcf.gz"
    mv results/variants/candidateSV.vcf.gz $eid_nameroot.candidateSV.vcf.gz
    all_candidates=$(dx upload "$eid_nameroot.candidateSV.vcf.gz" --brief)
    dx-jobutil-add-output all_candidates --class=file "$all_candidates"
    
    #File all_candidates_tbi = "results/variants/candidateSV.vcf.gz.tbi"
    mv results/variants/candidateSV.vcf.gz.tbi $eid_nameroot.candidateSV.vcf.gz.tbi
    all_candidates_tbi=$(dx upload "$eid_nameroot.candidateSV.vcf.gz.tbi" --brief)
    dx-jobutil-add-output all_candidates_tbi --class=file "$all_candidates_tbi"
    
    #File small_candidates = "results/variants/candidateSmallIndels.vcf.gz"
    mv results/variants/candidateSmallIndels.vcf.gz $eid_nameroot.candidateSmallIndels.vcf.gz
    small_candidates=$(dx upload "$eid_nameroot.candidateSmallIndels.vcf.gz" --brief)
    dx-jobutil-add-output small_candidates --class=file "$small_candidates"
    
    #File small_candidates_tbi = "results/variants/candidateSmallIndels.vcf.gz.tbi"
    mv results/variants/candidateSmallIndels.vcf.gz.tbi $eid_nameroot.candidateSmallIndels.vcf.gz.tbi
    small_candidates_tbi=$(dx upload "$eid_nameroot.candidateSmallIndels.vcf.gz.tbi" --brief)
    dx-jobutil-add-output small_candidates_tbi --class=file "$small_candidates_tbi"
    
    #File? tumor_only_variants = "results/variants/tumorSV.vcf.gz"
    mv results/variants/tumorSV.vcf.gz $eid_nameroot.tumorSV.vcf.gz
    tumor_only_variants=$(dx upload "$eid_nameroot.tumorSV.vcf.gz" --brief)
    dx-jobutil-add-output tumor_only_variants --class=file "$tumor_only_variants"
    
    #File? tumor_only_variants_tbi = "results/variants/tumorSV.vcf.gz.tbi"
    mv results/variants/tumorSV.vcf.gz.tbi $eid_nameroot.tumorSV.vcf.gz.tbi
    tumor_only_variants_tbi=$(dx upload "$eid_nameroot.tumorSV.vcf.gz.tbi" --brief)
    dx-jobutil-add-output tumor_only_variants_tbi --class=file "$tumor_only_variants_tbi"
    

}
