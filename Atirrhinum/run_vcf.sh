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
