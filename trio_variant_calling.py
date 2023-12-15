## Variant calling v2.0
import glob
import sys
import os
import re

############# Functions ##############
def get_data():
    arguments = sys.argv
    print('Arguments: %s' %arguments)
    data_folders = glob.glob('%s/%s' %(data_dir,arguments[1]))
    data_folders = [element for element in data_folders if element not in ('%s,%s')%(data_trimmed_dir,fastqc_dir)]
    print ('data_folders: %s,%s' %(len(data_folders),data_folders))
    return data_folders


def create_dirs(samtools_results,gatk_results,varscan_results,alignment_results,combined_variants):
    dirs = [samtools_results,gatk_results,varscan_results,alignment_results,combined_variants]
    for dir in dirs:
        # create results folder
        print(dir)
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
            print ('\033[31m %s directory doesn NOT exists. I am creating it. \033[0m' %(dir))
        else:
            print ('\033[31m %s directory exists. Not creating. \033[0m' %(dir))


################# samtools fixmate, sort and index to cleanup read pair info and flags #########################
def run_samtools_fixmate(alignment_results,bam_file,sample_id,files_2_delete):
    print
    print( "\033[34m Running SAMtools fixmate... \033[0m")
    
    # Create file names
    fixmate_file = os.path.join(alignment_results,sample_id + '_fixmate.bam')
    sorted_bam = os.path.join(alignment_results,sample_id + '_sorted.bam')
    
    # Define commands
    cmd1 = 'samtools fixmate -O bam %s %s' %(bam_file,fixmate_file)
    cmd2 = 'samtools sort  -@ 8 -O bam -o %s -T /proj/omics4tb/sturkarslan/tmp/%s_temp %s' %(sorted_bam, sample_id, fixmate_file)
    cmd3 = 'samtools index %s' %(sorted_bam)
    
    print("\n")
    print ("++++++ Samtools Fixmate Command: ", cmd1)
    #os.system(cmd1)
    print("\n")
    print ("++++++ Samtools Sort Command: ", cmd2)
    #os.system(cmd2)
    print("\n")
    print ("++++++ Samtools Index Command: ", cmd3)
    #os.system(cmd3)
    
    # add  temp files to list to delete
    temp_files = ['%s'%sorted_bam, '%s.bai'%sorted_bam, '%s'%fixmate_file]
    for temp_file in temp_files:
        files_2_delete.append(temp_file)


####################### Samtools Variant Calling ###############################
def samtools_variants(alignment_results,bam_file,samtools_results,sample_id): # With samtools
    print
    print( "\033[34m Running SAMtools Variant Calling.. \033[0m")
    
    # create file names
    sorted_bam = os.path.join(alignment_results,sample_id + '_sorted.bam')

    # Define commands
    cmd22 = 'samtools sort  -@ 8 -O bam -o %s -T /proj/omics4tb/sturkarslan/tmp/%s_temp %s' %(sorted_bam, sample_id, bam_file)
    cmd33 = 'samtools index %s' %(sorted_bam)
    
    print("\n")
    print ("++++++ Samtools Sort Command: \n", cmd22)
    os.system(cmd22)
    print("\n")
    print ("++++++ Samtools Index Command: \n", cmd33)
    os.system(cmd33)

    # create samtools results specific results directory
    samtools_files_path = '%s/%s'%(samtools_results,sample_id)

    # Produce BCF file with all locations in the genome
    cmd1 = 'bcftools  mpileup --threads 8 -x -O b -I -Q 13 -a "INFO/SCR,FORMAT/SP,INFO/ADR,INFO/ADF" -h 100 -M 10000 -f %s %s | bcftools call --ploidy 1 -O v -V indels -m -o %s_samtools.vcf' %(genome_fasta,sorted_bam, samtools_files_path)
    
    # Prepare vcf file for querying
    cmd2 = 'tabix -p vcf %s_samtools.vcf' %(samtools_files_path)
    
    # Variant Filtering
    percentageString = "%"
    cmd3 = " ~/bcftools-1.1/bin/bcftools filter -O z -o %s_samtools_final.vcf.gz -s LOWQUAL -i '%sQUAL>20 && DP>75 && VDB>0.00001 && MQ>30' %s_samtools.vcf" %(samtools_files_path, percentageString, samtools_files_path)
    
    print("\n")
    print( "++++++ Variant Calling mpileup: \n", cmd1)
    #os.system(cmd1)
    print
    print( "++++++ Variant Calling tabix: \n", cmd2)
    #os.system(cmd2)
    print
    print( "++++++ Variant Calling filtering: \n", cmd3)
    #os.system(cmd3)
    print
    
    #files_2_delete.append('%s_sorted.bam, %s_sorted.bam.bai, %s_fixmate.bam'%(base_file_name,base_file_name,base_file_name))
    return samtools_files_path,sorted_bam


####################### Varscan Variant Calling ###############################
def varscan_variants(sorted_bam,varscan_results,sample_id,files_2_delete): # with varscan
    print
    print( "\033[34m Running Varscan.. \033[0m")
    
    # create varscan results specific results directory
    varscan_files_path = '%s/%s'%(varscan_results,sample_id)

    # varscan mpileup
    cmd1 = 'samtools mpileup --input-fmt-option nthreads=8 -B -f %s -o %s.pileup %s' %(genome_fasta, varscan_files_path, sorted_bam)
    
    # varscan for snps
    cmd2 = '%s -Xmx128m -jar %s mpileup2snp %s.pileup --output-vcf 1 --min-coverage 8 --min-reads2 2 --min-avg-qual 30 --strand-filter 0 > %s_varscan_snps_final.vcf' %(javaPath, varscanPath, varscan_files_path, varscan_files_path)
    
    # varscan for indels
    cmd3 = '%s -Xmx128m -jar %s mpileup2indel %s.pileup --output-vcf 1 --min-coverage 8 --min-reads2 2 --min-avg-qual 30 --strand-filter 0 > %s_varscan_inds_final.vcf' %(javaPath, varscanPath, varscan_files_path, varscan_files_path)
    
    print("\n")
    print( "++++++ samtools Mpileup: \n", cmd1)
    print
    #os.system(cmd1)
    print("\n")
    print( "++++++ Varscan for SNPs: \n", cmd2)
    print
    #os.system(cmd2)
    print("\n")
    print
    print( "++++++ Varscan for INDELS: \n", cmd3)
    #os.system(cmd3)

    files_2_delete.append('%s.pileup' %(varscan_files_path))

    return varscan_files_path



####################### GATK Variant Calling ###############################
def gatk_variants(sorted_bam,gatk_results,sample_id): #with GATK HaploTypeCaller
    print
    print( "\033[34m Running GATK Haplotype Variant Caller.. \033[0m")
    
    # create varscan results specific results directory
    gatk_files_path = '%s/%s'%(gatk_results,sample_id)

    # haplotype command
    cmd1 = '%s -Xmx4G -jar %s -T HaplotypeCaller -R %s -I %s --genotyping_mode DISCOVERY -ploidy 1 -stand_emit_conf 30 -stand_call_conf 30 -o %s_gatk_raw.vcf' %(javaPath, gatkPath, genome_fasta, sorted_bam, gatk_files_path)
    
    # Select snp variants
    cmd2 = '%s -Xmx128m -jar %s -T SelectVariants -R %s -V %s_gatk_raw.vcf -selectType SNP -o %s_gatk_snps.vcf' %(javaPath, gatkPath, genome_fasta, gatk_files_path,gatk_files_path)
    
    # Apply filters to SNPs
    cmd3 = "%s -Xmx4G -jar %s -T VariantFiltration -R %s -V %s_gatk_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'my_snp_filter' -o %s_gatk_snps_filtered.vcf" %(javaPath, gatkPath, genome_fasta, gatk_files_path, gatk_files_path)
    
    # Select indel variants
    cmd4 = '%s -Xmx4G -jar %s -T SelectVariants -R %s -V %s_gatk_raw.vcf -selectType INDEL -o %s_gatk_inds.vcf' %(javaPath, gatkPath, genome_fasta, gatk_files_path,gatk_files_path)
    
    # Apply filters to indels
    cmd5 = "%s -Xmx4G -jar %s -T VariantFiltration -R %s -V %s_gatk_inds.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filterName 'my_indel_filter' -o %s_gatk_inds_filtered.vcf" %(javaPath, gatkPath, genome_fasta, gatk_files_path, gatk_files_path)
    
    # Merge vcf files
    cmd6 = "%s -Xmx4G -jar %s -T CombineVariants -R %s --variant %s_gatk_snps_filtered.vcf --variant %s_gatk_inds_filtered.vcf -o %s_gatk_final.vcf  -genotypeMergeOptions UNSORTED" %(javaPath, gatkPath, genome_fasta, gatk_files_path, gatk_files_path, gatk_files_path)

    print("\n")
    print( "++++++ GATK HaplotypeCaller Comnand: \n", cmd1)
    print
    os.system(cmd1)
    print("\n")
    print( "++++++ Select SNP Variants: \n", cmd2)
    print
    os.system(cmd2)
    print("\n")
    print( "++++++ Applying filters for SNPs: \n", cmd3)
    print
    os.system(cmd3)
    print("\n")
    print( "++++++ Select IndelVariants: \n", cmd4)
    print
    os.system(cmd4)
    print("\n")
    print( "++++++ Applying filters for Indelss: \n", cmd5)
    print
    os.system(cmd5)
    print("\n")
    print( "++++++ Merging vcf files: \n", cmd6)
    print
    os.system(cmd6)

    return gatk_files_path


####################### Run SNPEff annotations ###############################
def run_snpeff(samtools_files_path,varscan_files_path,gatk_files_path,combined_variants,sample_id):
    print
    print( "\033[34m Running SNPEff Annotations.. \033[0m")

    # Get list of VCF files from all callers
    vcf_files = ['%s_samtools_final.vcf.gz'%(samtools_files_path), '%s_varscan_inds_final.vcf'%(varscan_files_path), '%s_varscan_snps_final.vcf'%(varscan_files_path), '%s_gatk_final.vcf'%(gatk_files_path)]
    
    varscan_files = ['%s_varscan_inds_final.vcf'%(varscan_files_path), '%s_varscan_snps_final.vcf'%(varscan_files_path)]

    # create output file for combined variants output
    combined_variants_output = '%s/%s_combined_variants.txt' %(combined_variants,sample_id)
    print('combined_variants_output:%s' %(combined_variants_output))

    # open the final output file for writing before combining
    run_snpeff.t = open(combined_variants_output, 'w')

    # Loop through each VCF file
    for vcf_file in vcf_files:
        print(f'vcf_file{vcf_file}')
        if vcf_file == '%s_samtools_final.vcf.gz'%(samtools_files_path):
            # creat file names for output
            snpeff_vcf = re.split('final.', vcf_file)[0] + 'snpeff.vcf' 
            snpeff_filtered_vcf = re.split('final.', vcf_file)[0] + 'snpeff_filtered.vcf'
        else:
            # creat file names for output
            snpeff_vcf = re.split('final.', vcf_file)[0] + 'snpeff.' + re.split('final.', vcf_file)[1]
            snpeff_filtered_vcf = re.split('final.', vcf_file)[0] + 'snpeff_filtered.' + re.split('final.', vcf_file)[1]

        snpeff_stats = re.split('final.', vcf_file)[0] + 'snpeff_stats.txt'
        snpeff_final = re.split('final.', vcf_file)[0] + 'snpeff_final.txt'
        

        # snpeff formateff commands
        cmd1 ='%s -Xmx2g -jar %ssnpEff.jar -ud 0 -classic -csvStats %s -geneId -lof -v -formatEff -o gatk %s %s > %s'  %(javaPath, snpeff_path, snpeff_stats, snpeff_db, vcf_file, snpeff_vcf)
        
        # snpeff filtering command
        cmd2 = 'cat %s | %s -Xmx128m -jar %sSnpSift.jar filter "(FILTER = \'PASS\') & (EFF[*].CODING != \'NON_CODING\')" > %s' %(snpeff_vcf, javaPath, snpeff_path, snpeff_filtered_vcf)
       
        # create final one line variant file
        cmd3 ='cat %s | perl %sscripts/vcfEffOnePerLine.pl | %s -Xmx128m -jar %sSnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "(FILTER = \'PASS\')" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].CODING" "EFF[*].RANK" "EFF[*].DISTANCE"> %s'%(snpeff_filtered_vcf, snpeff_path, javaPath, snpeff_path, snpeff_final)

        print("\n")
        print( "++++++ Running SNPEff Formateff command: \n", cmd1)
        print
        os.system(cmd1)
        print("\n")
        print( "++++++ Running SNPEff Filtering: \n", cmd2)
        print
        os.system(cmd2)
        print("\n")
        print( "++++++ Running SNPEff Oneline final formatter: \n", cmd3)
        print
        os.system(cmd3)

        # Run function to combine vcf files into a single file from 3 callers
        combine_variants(snpeff_filtered_vcf,snpeff_final,combined_variants)
    run_snpeff.t.close()

####################### Combine variants from snpeff outputs ###############################
def is_file_empty(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as infile:
            line = infile.readline()
            return len(line) == 0
    else:
        return os.stat(path).st_size == 0

def vcfopen(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

####################### Combine variants from snpeff outputs ###############################
def combine_variants(snpeff_filtered_vcf, snpeff_final, combined_variants):
    """
    snpeff_filtered_vcf: input VCF file
    """
    print( "\033[34m Running Combine variants.. \033[0m")

    vcf_file_program = snpeff_filtered_vcf.split('/')[-2]
    print('snpeff_filtered_vcf:%s' %snpeff_filtered_vcf)

    ## open vcf file for processing and skip comment lines
    if is_file_empty(snpeff_filtered_vcf):
        print("FILE IS EMPTY !!!")
        return

    with vcfopen(snpeff_filtered_vcf) as infile:
        # skip all the comments
        prev_line = ''
        for line in infile:
            if not line.startswith("#"):
                break
            prev_line = line
        if not prev_line.startswith('#CHROM'):
            print("No chromosome information !!!, exiting combine_variants()")
            return

        frequencies = []
        for line in infile:
            vector = line.split("\t")

            ## if vcf file is varscan file
            if vcf_file_program == "varscan_results":
                ## collect vcf field for frequency
                freq = vector[9].split(':')[6]

                #print('++++++++Freq: %s+++++++++' %(freq))
                chm = line.split('\t')[0]
                my_fields = vector[7].split(';')
                #adp = myFields.split("ADP=")[1]
                ## append freq to list of freq
                frequencies.append(freq)

                ## if vcf file is samtools file
            elif vcf_file_program == "samtools_results":
                #collect list of all fields from vcf file
                my_fields = vector[7].split(';')
                #print(my_fields)
                # grab DP4 field
                dp4 = next(filter(lambda x:'DP4' in x, my_fields))
                #print dp4
                dp4_fields = dp4.split("DP4=")[1]
                fr = float(dp4_fields.split(",")[0])
                rr = float(dp4_fields.split(",")[1])
                fa = float(dp4_fields.split(",")[2])
                ra = float(dp4_fields.split(",")[3])
                freq = (fa + ra) / (fr + rr + fa + ra)
                freq = round(freq * 100, 2)
                freq = str(freq) + "%"

                frequencies.append(freq)

                ## if vcf file is gatk file
            elif vcf_file_program == "gatk_results":
                ## collect vcf field
                my_fields = vector[7].split(':')
                #print myFields
                if my_fields[0] != ".":
                    ad = float(my_fields[1].split(',')[1])
                    dp = float(my_fields[2])
                    freq = ad / dp
                    freq = round(freq * 100, 2)
                    freq = str(freq) + "%"
                ## append freq to list of freq
                    frequencies.append(freq)
                else:
                    frequencies.append('')
            else:
                frequencies.append('')

        print( 'Length of Frequencies= '+ str(len(frequencies)))

        # filename for the output from converting vcf to oneliner with frequency and program added
        outfile_w_freq = combined_variants + '/' + snpeff_filtered_vcf.split('/')[-1].split('_snpeff_filtered')[0] + '_outfile_w_freq.txt'

        with open(outfile_w_freq, 'w') as outfile, open(snpeff_final, 'r') as final_in:
            print(snpeff_final)
            final_in.readline()  # skip header

            myList = []
            index = 0
            for line in final_in:
                
                # change alternative chromosome names
                chromosome = line.split('\t')[0]
                if chromosome == "NC_002937":
                    print( "Found alternative chromosome name: %s" %chromosome)
                    chm = "Chromosome"
                elif chromosome == "NC_005863":
                    print( "Found alternative chromosome name: %s" %chromosome)
                    chm = "pDV"
                else:
                    chm = line.split('\t')[0]
                
                # collect fields
                if index < len(frequencies):
                    print(f'index: {frequencies[index]}')
                    pos = line.split('\t')[1]
                    refs = list(line.split('\t')[2])[0] # grab only first character in reference sequence
                    ref = line.split('\t')[2]
                    alt = line.split('\t')[3]
                    eff = line.split('\t')[9]
                    imp = line.split('\t')[10]
                    fnc = line.split('\t')[11]
                    cdn = line.split('\t')[12]
                    aac = line.split('\t')[13]
                    loc = line.split('\t')[15]
                    cod = line.split('\t')[16]
                    dep = line.split('\t')[6]
                    fre = frequencies[index]
                    pro = vcf_file_program.split('_')[0]

                    lineToWrite = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chm,pos,refs,ref,alt,eff,imp,fnc,cdn,aac,loc,cod,fre,pro,dep)
                    outfile.write(lineToWrite)

                    myList.append(lineToWrite)
                    index += 1
                else:
                    print("Index is out of range.")

    # write into final output file
    for element in myList:
        run_snpeff.t.write('%s' %element)

####################### Collate variants from three programs ###############################
def salomon(aggregatedList):
    #this function takes a list of lists of variants and find the consensus list

    # f.1. finding the unique positions
    uniqueLocations=[]
    for variants in aggregatedList:
        for variant in variants:
            uniqueLocation=variant[:3]
            if uniqueLocation not in uniqueLocations:
                uniqueLocations.append(uniqueLocation)

    # f.2. building the full consensus list
    consensus_list=[]
    for uniqueLocation in uniqueLocations:
        callers=[]
        body=[]
        freqs = []
        dps = []
        freq=''
        dp =''
        freqFloats = []

        for variants in aggregatedList:
            for variant in variants:
                if uniqueLocation == variant[:3]:
                    body=variant[:-2]
                    callers.append(variant[-2])
                    if variant[-2] == 'varscan':
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)
                        #print variant[:-2]
                        dp = "NA"
                        dps.append(dp)

                        if freq == '':
                            print( 'WARNING varscan did not provide frequency value for variant')
                            stringVariant='\t'.join(variant)
                            #print stringVariant

                    if variant[-2] == 'samtools':
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)
                        dp =variant[-1]
                        dps.append(dp)
                        if freq == '':
                            print( 'WARNING samtools did not provide frequency value for variant')
                            stringVariant='\t'.join(variant)
                            #print stringVariant

                    if variant[-2] == 'gatk':
                        dp = variant[-1]
                        dps.append(dp)
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)

        # incorporating what we found
        callersString=':'.join(callers)
        freqsString=':'.join(freqs)
        dpsString=':'.join(dps)
        freqFloatString = max(freqFloats)

        # Get the samtools frequency as final frequency if available otherwise pick varscan freq.
        if 'samtools' in callers and len(callers) > 1:
            final_freq = freqs[callers.index('samtools')]
        elif 'varscan' in callers and 'samtools' not in callers:
            final_freq = freqs[callers.index('varscan')]
        else:
            final_freq = freqs[0]

        body.append(callersString)
        body.append(freqsString)
        body.append(dpsString)
        body.append(freqFloatString)
        body.append(final_freq)
        consensus_list.append(body)

    return consensus_list


####################### Variant Retriever ###############################
def variantRetriever(combined_output_file):
    print(f'combined_output_file: {combined_output_file}')

    # this function retrieves the variants for each caller
    varscan_list = []
    gatk_list = []
    samtools_list = []
    print('combined_output')
    print(os.stat(combined_output_file).st_size)

    # Start reading each line and append to appropriate program list
    with open(combined_output_file) as f:
        for line in f:
            vector = line.split('\t')
            vector[-1]=vector[-1].replace('\n','')
            program=vector[-2]
            if program == 'varscan':
                varscan_list.append(vector)
            if program == 'gatk':
                gatk_list.append(vector)
            if program == 'samtools':
                samtools_list.append(vector)
        f.close()

    return varscan_list,gatk_list,samtools_list


####################### Collate variants ###############################
def collate_variants(combined_output_file,merged_variants_file):
    print
    print( "\033[34m Running Collate variants.. \033[0m")

    # 2. recovering list of variants
    varscan_list,gatk_list,samtools_list = variantRetriever(combined_output_file)
    print( 'detected variants',len(varscan_list),len(gatk_list),len(samtools_list))

    # 3. finding consensus list of variants
    print( 'merging...')
    consensus_list = salomon([varscan_list,gatk_list,samtools_list])
    print( 'final set',len(consensus_list))

    # 4. writing a consensus list
    print( 'writing file...')
    #print(merged_variants_file)

    g = open(merged_variants_file, 'w')
    for element in consensus_list:
        line2write='\t'.join(element)
        line2write=line2write+'\n'
        g.write(line2write)
    g.close()
    #return collate_variants.combined_output_file


####################### Delete temporary files ###############################
def delete_temp_files(files_2_delete):
    print
    print( "\033[34m Deleting Temporry files.. \033[0m")
    for file in files_2_delete:
        cmd = 'rm %s' %(file)
        print(cmd)
        #os.system(cmd)


####################### Running the Pipeline ###############################
def run_pipeline(organism,project_dir,experiment_name,sample_id):
    files_2_delete = []
    
    # Find the bam file
    results_dir = os.path.join(project_dir,"output",experiment_name,organism,sample_id)
    basecalls_dir = os.path.join(project_dir,"output",experiment_name,"basecalls",sample_id)
    path = os.path.join(basecalls_dir,"**/*.bam")
    bam_file = glob.glob(path, recursive=True)[0]
    
    # Creating program-specific results directories
    samtools_results = os.path.join(results_dir, 'samtools_results')
    gatk_results = os.path.join(results_dir, 'gatk_results')
    varscan_results = os.path.join(results_dir, 'varscan_results')
    alignment_results = os.path.join(results_dir, 'alignment_results')
    combined_variants = os.path.join(results_dir, 'combined_variants')

    # Creating final results files
    combined_output_file = os.path.join(combined_variants, f'{sample_id}_combined_variants.txt')
    merged_variants_file = os.path.join(combined_variants, f'{sample_id}_merged_variants_final.txt')

    base_file_name = os.path.join(alignment_results,sample_id)

    # Print paths for verification
    print(f'Organism: {organism}')
    print(f'sample_name: {sample_id}')
    print(f'experiment_name: {experiment_name}')
    print(f'results_dir: {results_dir}')
    print(f'BAM File: {bam_file}')
    print(f'samtools_results: {samtools_results}')
    print(f'gatk_results: {gatk_results}')
    print(f'varscan_results: {varscan_results}')
    print(f'alignment_results: {alignment_results}')
    print(f'combined_variants: {combined_variants}')
    print(f'combined_output_file: {combined_output_file}')
    print(f'merged_variants_file: {merged_variants_file}')

    
    # 01. Get directories
    #create_dirs(samtools_results,gatk_results,varscan_results,alignment_results,combined_variants)
    
    # 02. Run samtools fixmate
    #run_samtools_fixmate(alignment_results,bam_file,sample_id,files_2_delete)
    
    # 03. Run Samtools variant calling
    samtools_files_path,sorted_bam = samtools_variants(alignment_results,bam_file,samtools_results,sample_id)
    
    # 04. Run varscan variant calling
    varscan_files_path = varscan_variants(sorted_bam,varscan_results,sample_id,files_2_delete)
    
    # 05. Run GATK variant Calling
    gatk_files_path = gatk_variants(sorted_bam,gatk_results,sample_id)
    
    # 06. Run SNPEff annotations
    vcf_file = run_snpeff(samtools_files_path,varscan_files_path,gatk_files_path,combined_variants,sample_id)
    
    # 07. Collate variants into single file from 3 callers and unify them
    collate_variants(combined_output_file,merged_variants_file)
    
    # 08. Delete temporary files_2_delete
    #delete_temp_files(files_2_delete)

############# Programs #############
javaPath = "/usr/bin/java"
piccardPath = "/users/sturkars/picard-tools-1.139/picard.jar"
gatkPath = "/users/sturkars/gatk/GenomeAnalysisTK.jar"
varscanPath = "/users/sturkars/VarScan.v2.3.9.jar"
snpeff_path = "/users/sturkars/snpEffv43/snpEff/"

############# Globals ##############
organism = "dvh"
experiment_name = "EPD9"
sample_id = 'barcode07'

# data and results directories
project_dir = "/proj/omics4tb2/sturkarslan/Nanopore"
genome_dir = os.path.join(project_dir,"reference")
genomeGff = os.path.join(genome_dir,organism + ".GCA_000195755.1.30.gtf")
known_sites = os.path.join(project_dir,"reference", organism + "-variants-compiled_sorted.vcf")

######### Annotation databases #######
# snpEff databases
if organism == "mmp":
    snpeff_db = "Methanococcus_maripaludis_s2"
    genome_fasta = os.path.join(genome_dir,'DvH_Mmp_merged_dna.genome.fasta')
    #genome_fasta = os.path.join(genome_dir,'Methanococcus_maripaludis_s2.GCA_000011585.1.30.dna.genome.fasta')
if organism == "dvh":
    snpeff_db = "Desulfovibrio_vulgaris_str_hildenborough"
    genome_fasta = os.path.join(genome_dir,'DvH_Mmp_merged_dna.genome.fasta')
    #genome_fasta = os.path.join(genome_dir,'Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta')


run_pipeline(organism, project_dir, experiment_name, sample_id)
