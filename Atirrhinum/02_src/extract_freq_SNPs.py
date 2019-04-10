#!/usr/local/bin/python

from __future__ import division
import sys
import gzip
import copy

vcf_file = gzip.open(sys.argv[1])
name_file = open(sys.argv[2])
order_sample = open('01_raw_data/vcf_file_order.list')

name_dict = {}
name_list = []

for i in name_file:
	i = i.strip().split()
	if i[2] not in name_dict:
		name_dict[i[2]] = [int(i[0].split('_')[1])]
		name_list.append(int(i[0].split('_')[1]))
	else:
		name_dict[i[2]].append(int(i[0].split('_')[1]))
		name_list.append(int(i[0].split('_')[1]))
### name_dict = {A.hispanicum_Morocco [102, 59, 58], A.graniticum [25, 26, 27, 28, 29, 30, 101], A.boissieri [7, 8, 9]...}

position_in_vcf = []
for i in order_sample:
	i = i.strip().split()
	for a in i:
		if '.bam' in a:
			a = a.split('.')
			position_in_vcf.append(int(a[0][7:]))
		else:
			position_in_vcf.append('')
## position_in_vcf = ['', '', '', '', '', '', '', '', '', 100, 101, 102, 104, 105, 106, 107, 108, 109, 10, 110...]
def allele_freq(species_alleles):
	freq = {0:0,1:0}
	for ind_allele in species_alleles:
		if ind_allele == '':
			continue
		else:
			if len(ind_allele) == 1:
				if ind_allele == '0':
					freq[0] += 2 ## homozygous in first allele
				else:
					freq[1] += 2
			else:
				freq[0] += 1
				freq[1] += 1
	return freq

#out_file_Pairwise_Allele_Freq = open('Pairwise_Allele_Freq.mid.txt','w')
#out_file_Pairwise_Allele_Freq.write('Pop_A_name\tPop_A_Allele_Freq\tPop_B_name\tPop_B_Allele_Freq\tAll_Allele_Freq(exp_A)\n')

# for every SNP locus, put all individuals's allele in a list: [0,1,1,0,01,12,1,12,2,012....]
for i in vcf_file:
	if not i.startswith('#'):
		all_alleles = ['']*120
		i = i.strip().split()
		for order in range(9,110):
			sample = i[order].split(':')
			alleles = sample[-1].split(',')
			for alt in range(0,len(alleles)):
				if int(alleles[alt]) > 0:
					all_alleles[order] = all_alleles[order]+str(alt)
# for every species, put its individuals into a list, the rest another list, create a list for every other species and put them into a dictionary-{rest_species_allele_dict}
		for A in name_dict:
			rest_alleles = []
			species_alleles = []
			rest_individual_list, rest_species_name, rest_species_allele_dict = name_list[:], copy.deepcopy(name_dict),{}
			del rest_species_name[A]
			FREQ_DIFF_out_file = open('./03_output/'+A+'.FREQ_DIFF.list','a')
			for individuals in name_dict[A]:
				rest_individual_list.remove(individuals)
				species_alleles.append(all_alleles[position_in_vcf.index(individuals)])
			for rest_individuals in rest_individual_list:
				rest_alleles.append(all_alleles[position_in_vcf.index(rest_individuals)])
			freq_A, freq_ALL = allele_freq(species_alleles), allele_freq(rest_alleles)
			rest_group_freq_0_list,rest_group_freq_1_list = [],[]
			for B in rest_species_name:
				B_species_alleles = []
				for B_individuals in rest_species_name[B]:
					B_species_alleles.append(all_alleles[position_in_vcf.index(B_individuals)])
				rest_species_allele_dict[B] = B_species_alleles
				freq_B = allele_freq(B_species_alleles)
				if freq_B[0]+freq_B[1] != 0:
					AF_B_0 = freq_B[0]/(freq_B[0]+freq_B[1])
					AF_B_1 = freq_B[1]/(freq_B[0]+freq_B[1])
					rest_group_freq_0_list.append(AF_B_0)
					rest_group_freq_1_list.append(AF_B_1)
# judge if alleles_frequency in species A is significantly different from all the rest species respetively and in all.
## retaining loci by thresholds
			if len(rest_group_freq_0_list) > 5: ## for this locus, alleles in at least 5 other groups are not missing
				if freq_A[0]+freq_A[1] != 0 and freq_ALL[0]+freq_ALL[1] != 0:
					AF_A_0 = freq_A[0]/(freq_A[0]+freq_A[1])
					AF_A_1 = freq_A[1]/(freq_A[0]+freq_A[1])
					AF_ALL_0 = freq_ALL[0]/(freq_ALL[0]+freq_ALL[1])
					AF_ALL_1 = freq_ALL[1]/(freq_ALL[0]+freq_ALL[1])
# stage 1:  If AF in this group is higher than a specific threshold (say 90%), progresses to stage 2 ->
					if AF_A_0 > 0.9:
# stage 2: if AF in all the rest groups(aggregated) is lower than a threshold (10%), progress to stage 3 ->
						if AF_ALL_0 <= 0.1:
# stage 3: if AFs in any other groups (separated) are lower than a threshold (10%), take this locus as valid one with allele frequency difference
							if AF_A_0 == 1 and AF_ALL_0 == 0:
								continue
							else:
								if max(rest_group_freq_0_list) < 0.1:
									FREQ_DIFF_out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(i[0], i[1], freq_A, round(AF_A_0,2), round(AF_A_1,2), round(AF_ALL_0,2), round(AF_ALL_1,2), round(max(rest_group_freq_0_list),2)))
					elif AF_A_1 > 0.9:
# stage 2: if AF in all the rest groups(aggregated) is lower than a threshold (10%), progress to stage 3 ->
						if AF_ALL_1 <= 0.1:
							if AF_A_1 == 1 and AF_ALL_1 == 0:
								continue
							else:
								if max(rest_group_freq_0_list) < 0.1:
									FREQ_DIFF_out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(i[0], i[1], freq_A, round(AF_A_0,2), round(AF_A_1,2), round(AF_ALL_0,2), round(AF_ALL_1,2), round(max(rest_group_freq_0_list),2)))
#								print '03',AF_A_0, AF_A_1, A,freq_A, freq_ALL, max(rest_group_freq_0_list), rest_group_freq_0_list
					else:
						continue
#					out_file_Pairwise_Allele_Freq.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(A,freq_A, B, freq_B, freq_ALL))

vcf_file.close()
name_file.close()
