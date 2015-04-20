import sys, gzip 
from JTfunctions import * 

dose_file = sys.argv[1]
blup_file = sys.argv[2]
out_file = sys.argv[3]

snplist = list_from_file(blup_file,1)
snpdic = {}
fin = open(blup_file,'r')
for line in fin:
	l = line.strip().split()
	snp, eff = l[0],l[2]
	snpdic[snp] = eff
fin.close()

fout = gzip.open(out_file,'wb')

din = gzip.open(dose_file, 'rb')
for line in din:
	l = line.strip().split()
	names = l.pop(0)
	write_list = names.split("->")
	l.pop(0)
	sum = 0 
	if len(l) == len(snplist):
		for i in range(0,len(l)):
			dosage = float(l[i])
			effect = float(snpdic[snplist[i]])
			score = dosage * effect
			sum += score 
	else:
		print "FAIL"
	write_list = write_list + [str(sum)]
	fout.write("\t".join(write_list)+"\n")
din.close()
fout.close()





