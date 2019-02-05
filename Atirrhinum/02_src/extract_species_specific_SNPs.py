#!/usr/local/bin/python

import sys
import gzip

vcf_file = gzip.open(sys.argv[1])
name_file = open(sys.argv[2])
out_file = open(sys.argv[3],'w')
order_sample = open('../01_raw_data/vcf_file_order.list')

name_dict = {}

for i in name_file:
	i = i.strip().split()
	if i[2] not in name_dict:
		name_dict[i[2]] = [int(i[0].split('_')[1])]
	else:
		name_dict[i[2]].append(int(i[0].split('_')[1]))

position_in_vcf = []
for i in order_sample:
	i = i.strip().split()
	for a in i:
		if '.bam' in a:
			a = a.split('.')
			position_in_vcf.append(int(a[0][7:]))
		else:
			position_in_vcf.append('')

#print position_in_vcf

#for i in name_dict:
#	print i,name_dict[i]

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
#		print all_alleles

		rest_alleles = []
		species_alleles = []
		for A in name_dict:
			if A == 'A.boissieri':
				for individuals in name_dict[A]:
					species_alleles.append(all_alleles[position_in_vcf.index(individuals)])
#				print 'the allele status of ',A,' is: ',species_alleles  ## the allele status of  A.boissieri  is:  ['1', '1', '1']
			else:
				for individuals in name_dict[A]:
					rest_alleles.append(all_alleles[position_in_vcf.index(individuals)])
#		print  'the allele status of the rest is: ', rest_alleles

		if len(set(species_alleles)) == 1 and set(species_alleles) != set(species_alleles):
			if set(rest_alleles) > set(species_alleles):  ## alleles shows not only in A
				continue
			else:
				print 'the allele status of ',A,' is: ',species_alleles  ## the allele status of  A.boissieri  is:  ['1', '1', '1']
				print  'the allele status of the rest is: ', rest_alleles
				out_file.write('{0}\t{1}\t{2}\tSSL\n'.format(i[0],i[1],A))
'''
		else:
			print frequency
'''

vcf_file.close()
name_file.close()
out_file.close()
