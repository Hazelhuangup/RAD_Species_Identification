#!/usr/local/bin/python
import sys

f1 = open(sys.argv[1])
f2 = open(sys.argv[2])

d = {}
for i in f1:
	if i.startswith('sed'):
		i = i.split(' ')
		a = i[2].split('/')
		d[a[2]] = a[1][2:-2]

for i in f2:
	i = i.split()
	print d[i[0]],i[1]

f1.close
f2.close
