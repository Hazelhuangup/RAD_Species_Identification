#!/usr/local/bin/python

import sys
import gzip

vcf_file = gzip.open(sys.argv[1])
name_file = open(sys.argv[2])
order_sample = open('./01_raw_data/vcf_file_order.list')

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

position_in_vcf = []
for i in order_sample:
	i = i.strip().split()
	for a in i:
		if '.bam' in a:
			a = a.split('.')
			position_in_vcf.append(int(a[0][7:]))
		else:
			position_in_vcf.append('')

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
# for every species, put its individuals into a list, the rest another list.
		for A in name_dict:
			rest_alleles = []
			species_alleles = []
			rest_individual_list = name_list[:]
			out_file = open('./03_output/'+A+'.SSA.list','a')
			for individuals in name_dict[A]:
				rest_individual_list.remove(individuals)
				species_alleles.append(all_alleles[position_in_vcf.index(individuals)])
			for individuals in rest_individual_list:
				rest_alleles.append(all_alleles[position_in_vcf.index(individuals)])
# judge if all individuals in species A are consensus and different from all the rest species.
			if len(set(species_alleles)) == 1 and set(species_alleles) != set(rest_alleles):
				if set(rest_alleles) > set(species_alleles):  ## alleles shows not only in A
					continue
				else:
					out_file.write('{0}\t{1}\t{2}\tSSL\n'.format(i[0],i[1],A))

for A in name_dict:
	file_name = A+'.SSA.list'
	file_name.close()

'''
		else:
			print frequency
'''

vcf_file.close()
name_file.close()
