After the GENOMICON-Seq runs were finished, generated FASTQ files were mapped with BWA (`bwa_mapping.sh`), while produced BAM files are also sorted. 

Generated BAM files were processed with Strelka2 (`strelka2.sh`). The script compares the control sample without any mutations but having the same number of reads and the tumor samples, and tumor sample with the inserted mutation that mimics the chosen SBS-signature. `strelka2.slurm` run the `strelka2.sh` on HPC. 

