## Count the depth on the whole genome
vcftools --gzvcf ./01_raw_data/allsamples121217_raw_m.vcf.gz --freq --out vcf_basic_info
vcftools --gzvcf ./01_raw_data/allsamples121217_raw_m.vcf.gz --site-depth --out 03_output/vcf_site_depth
vcftools --gzvcf ./01_raw_data/allsamples121217_raw_m.vcf.gz --depth --out 03_output/vcf_indi_depth
vcftools --gzvcf ./01_raw_data/allsamples121217_raw_m.vcf.gz --geno-depth --out 03_output/vcf_indi_depth
vcftools --gzvcf ./01_raw_data/allsamples121217_raw_m.vcf.gz --SNPdensity 10000 --out 03_output/vcf_SNPdensity

## Only keep SNP/INDEL loci
zcat ./01_raw_data/allsamples121217_raw_m.vcf.gz |awk '{if(/^##/) print $0;else{if($5!=".") print $0}}' |gzip -c > ./01_raw_data/allsamples121217_var.vcf.gz
## Only keep SNP loci
vcftools --gzvcf ./01_raw_data/allsamples121217_var.vcf.gz --remove-indels --recode --recode-INFO-all --stout|gzip -c > ./01_raw_data/allsamples121217_SNPs_only.vcf.gz ## 37min
## Only keep Bi-allelic SNP loci
zcat ./01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz |awk '{if(/^##/) print $0;else{if(length($5)==1) print $0}}' |gzip -c > ./01_raw_data/allsamples121217_Bi-allelic_SNP_only.vcf.gz

## extract hetero-loci in sample 1
zcat ./01_raw_data/allsamples121217_var.SNPs_only.vcf.gz|awk '{split($10,a,":");split(a[7],b,",");if(b[1]!=0 && b[2]!=0) print $1"\t"$2"\t"$4"\t"$5"\t"$10}' > sample_100.hetero_loci

## re_correspond the species name to samples
./02_src/rename_samples.py ./01_raw_data/number2names.sh ./01_raw_data/antirrhinum_names.txt > ./01_raw_data/Sample_to_antirrhinum_names.txt

## testing the species specific alleles in A.bolssierri (sample 7/8/9)
./02_src/extract_species_specific_SNPs.py  ./01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz ./01_raw_data/Sample_to_antirrhinum_names.txt ./01_raw_data/A.bolssierri.SSA

## extract species specific alleles in all species
./02_src/extract_species_specific_SNPs.py  ./01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz ./01_raw_data/Sample_to_antirrhinum_names.txt 

## basic statistical analysis and visualization of SSA
for i in `ls ./03_output/A.*list`;do
            wc -l $i >> ./03_output/Atirrhinum.SSA_number_for_every_species.list
            if [ -s $i ]; then
                        ./02_src/sliding_window.py $i $i\.sliding_win.freq
            fi
done
sort -n ./03_output/Atirrhinum.SSA_number_for_every_species.list

Rscript ./02_src/SSA_frequency_sliding_windows_on_chr.R ./03_output/

### calculating allele frequency difference
./02_src/extract_freq_SNPs.py ./01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz ./01_raw_data/Sample_to_antirrhinum_names.txt 

for i in `ls 03_output/*FREQ_DIFF.list`;do
	wc -l $i;
	for z in {1..8};do
		grep 'Chr'$z $i|wc -l 
		done
	done > 03_output/All.FREQ_DIFF.No.txt

for i in {10..110};do
	zcat ./01_raw_data/allsamples121217_Bi-allelic_SNP_only.vcf.gz | awk '{if(!/^##/) split($'$i',a,":"); print a[3]}' |sort -n |uniq -c  > 03_output/read_depth_for_ind_freq/$i.freq
	done

Rscript ./02_src/draw_ind_depth_frequency.R /localdisk/home/s1847348/02_project/01_nuclear_for_barcoding

./02_src/count_missing_individuals_loci.py ./01_raw_data/allsamples121217_Bi-allelic_SNP_only.vcf.gz ./01_raw_data/allsamples121217_Bi-allelic_SNP_only.count_missing_ind.list
awk '{print $3}' ./01_raw_data/allsamples121217_Bi-allelic_SNP_only.count_missing_ind.list|sort -n |uniq -c |awk '{print $1"\t"$2}' > draw_frequency_no_ind_at_a_locus.freq
Rscript ./02_src/draw_frequency_basic.R numbers_of_ind_not_missing_at_this_locus.pdf draw_frequency_no_ind_at_a_locus.freq

### calculating Fst
### extracting information from vcf to -> .tfam, .tped
#vcf2tped.py 01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz allsamples121217_SNPs_only.vcf.recode.vcf
#
### rename the .tfam file using the real taxonomy name
#python 02_src/modify_names_in_tfam_file.py 03_output/allsamples121217_SNPs_only.vcf.recode.vcf.tfam 01_raw_data/Sample_to_antirrhinum_names.txt  03_output/allsamples121217_SNPs_only.vcf.recode.vcf.mod.tfam
#
### calculte Fst using popstats.py
#vcftools --gzvcf ./01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz --weir-fst-pop A.barrelieri.pop --weir-fst-pop A.graniticum.pop --out barrelieri_vs_graniticum
#vcftools --gzvcf ./01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz --weir-fst-pop A.barrelieri.pop --weir-fst-pop A.tortuosum --out barrelieri_vs_tortuosum
#vcftools --gzvcf ./01_raw_data/allsamples121217_SNPs_only.vcf.recode.vcf.gz --weir-fst-pop A.tortuosum.pop --weir-fst-pop A.graniticum.pop --out tortuosum_vs_graniticum
