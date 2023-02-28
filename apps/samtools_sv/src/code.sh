#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail

    dx-download-all-inputs --parallel

    echo ${cram_path}
    echo ${cram_index_path}
    echo ${fusion_sites_path}
    echo ${reference_path}
    echo ${reference_fai_path}
    echo ${dockerimage_samtools_sv_path}

    mv $cram_index_path ~/in/cram
    mv $reference_fai_path ~/in/reference


    docker load -i ${dockerimage_samtools_sv_path}
    ## get <eid>_23153_0_0 from <eid>_23153_0_0.bqsr.bam
    eid_nameroot=$(echo $cram_name | cut -d'.' -f1)
    
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /usr/local/reference/:/usr/local/reference/ -w /home/dnanexus chrisamiller/samtools_sv_filter:0.1 \
        /bin/bash /usr/local/bin/samtools_sv.sh ${cram_path} ${reference_path} ${fusion_sites_path} ${eid_nameroot} 

    # subset_bam="$eid_nameroot.subset.bam"
    subset_bam=$(dx upload "$eid_nameroot.subset.bam" --brief)
    dx-jobutil-add-output subset_bam --class=file "$subset_bam"
    
    # subset_bam_bai="$eid_nameroot.subset.bam.bai"
    subset_bam_bai=$(dx upload "$eid_nameroot.subset.bam.bai" --brief)
    dx-jobutil-add-output subset_bam_bai --class=file "$subset_bam_bai"
    
    # subset_filtered_bam="$eid_nameroot.subset_filtered.bam"
    subset_filtered_bam=$(dx upload "$eid_nameroot.subset_filtered.bam" --brief)
    dx-jobutil-add-output subset_filtered_bam --class=file "$subset_filtered_bam"
    
    # subset_filtered_bam_bai="$eid_nameroot.subset_filtered.bam.bai"
    subset_filtered_bam_bai=$(dx upload "$eid_nameroot.subset_filtered.bam.bai" --brief)
    dx-jobutil-add-output subset_filtered_bam_bai --class=file "$subset_filtered_bam_bai"
    
    # flagstat="$eid_nameroot.allreads.flagstat"
    flagstat=$(dx upload "$eid_nameroot.allreads.flagstat" --brief)
    dx-jobutil-add-output flagstat --class=file "$flagstat"
    
    # readlength="$eid_nameroot.read_lengths.txt"
    readlength=$(dx upload "$eid_nameroot.read_lengths.txt" --brief)
    dx-jobutil-add-output readlength --class=file "$readlength"
    

}
