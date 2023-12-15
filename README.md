# NanoVar-Micro
 Nanopore Microbial Variant Calling Pipeline

This pipeline runs dorado for high accuracy basecalling and alignment to reference genomes creating a BAM file that can be used as an input the bwa_pipeline_nanopore. This pipeline makes variant calls by using samtools/bcftools, varscan and gatk. It annotates the variants by using snpEff/snpSift and combines/collates variants that can be filtered based on the number of callers. 

# Setup of the environment
Basecalling with **dorado** requires GPU. Basecalling/alignment steps were performed in Macboo Pro (M1),3.2 Ghz, 10 CPU/16 GPU, 32GB Memory.

It took about ~13 hours to make basecalling and alignment to merged Dvh and Mmp genomes by using high quality basecalling model.

# Setup directory structure
- Nanopore
    - code
    - data
        - EPD9
            - barcode07
                - Run1, Run2, Run3, Run5
    - output
        - EPD9
            - dvh
                - barcode07
                    - alignment_results
                    - samtools_results
                    - gatk_results
                    - varscan_results
                    - combined_output
            - barcode07
                - basecalls
    - model
    - reference

# Install dorado
Download provided binaries for the relevant platform from [dorado Github repository](https://github.com/nanoporetech/dorado).

Unpack the dorado under the **code** directory

# Download the relevant dorado model
Download the relevant model from dorado gihub repository. We have used dna_r10.4.1_e8.2_400bps_sup@v4.3.0 model for the latest high accuracy basecalling.

    code/dorado-0.4.3-osx-arm64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0

model file is placed under the **model** directory

# Setting up POD5 for merging and converting fast5 files
    pip install pod5


# Convert .fats5 files from all the runs into one final .pod5 file
    pod5 convert fast5 -r Run1/barcode07/*.fast5 Run2/barcode07/*.fast5 Run3/barcode07/*.fast5 Run4/barcode07/*.fast5 --output barcode07_combined.pod5

# Run dorado
    code/dorado-0.4.3-osx-arm64/bin/dorado basecaller model/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 data/EPD9/barcode07/barcode07_combined.pod5 --reference reference/DvH_Mmp_merged_dna.genome.fasta > output/EPD9/barcode07/basecalls/barcode07.bam

# Running trio_variant_calling pipeline
Trio Variant Calling Pipeline uses BAM input file from dorado and calls variant by using samtools, gatk and varscan. Annotations are performed with snpEFF and variants are collated.
- Requirements
    - samtools
    - bcftools
    - varscan
    - gatk
    - snpEff
    
        
