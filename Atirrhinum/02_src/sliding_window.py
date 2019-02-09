#!/usr/bin/python2.7
import sys
SSA_file = open(sys.argv[1])
Sliding_window_freq_file = open(sys.argv[2],'w')

dict_SSA, n = {}, 0
for i in SSA_file:
    i = i.strip().split()
    if i[0] not in dict_SSA:
        dict_SSA[i[0]] = [int(i[1])]
    else:
        dict_SSA[i[0]].append(int(i[1]))

for i in dict_SSA:
    chr_length = max(dict_SSA[i])
    for window in range(0,int(chr_length/1000000)+2):
        start = window * 1000000
        end = start + 10000000 ### 10M step,1M SLIDING WINDOW
        for SSA in dict_SSA[i]:
            if SSA>=start and SSA<end:
                n += 1
        Sliding_window_freq_file.write('{0}\t{1}\t{2}\t{3}\n'.format(i,start,end,n))
        n = 0

SSA_file.close
Sliding_window_freq_file.close

