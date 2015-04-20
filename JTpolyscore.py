#!/python/bin/python -O 
# Jason Matthew Torres

'''
Python script to run Polyscore Procedure 
Requires: GCTA 1.24.4 
'''
#libraries 
import sys, os, gzip 
import subprocess as sp
from JTfunctions import * 
#globals 
cout = sys.stdout.write
cerr = sys.stderr.write

mach_dir = "/group/im-lab/nas40t2/jason/projects/PrediXmod/PrediXcan/StarrCounty/dosages/"
main_dir = "/group/im-lab/nas40t2/jason/projects/Polyscore_StarrCounty/"
# chr7.mldose.gz
grm_dir, pheno_dir, reml_dir, blup_dir = main_dir+"grms/", main_dir+"pheno/", main_dir+"reml/", main_dir+"blup/"
poly_dir = main_dir+"polyscore/"
#functions 

def make_grm(chrom):
	'''
	GCTA GRM function that returns name of grm 
	'''
	global mach_dir, grm_dir 
	dname, iname = "chr"+str(chrom)+".mldose.gz" , "chr"+str(chrom)+".mlinfo.gz" 
	dose, info  = mach_dir+dname, mach_dir+iname 
	dose_cp, info_cp = grm_dir+dname, grm_dir+iname
	cout("Copying dose and info files for chromosome: %s\n" % str(chrom))
	call1 = "cp %s %s" % (info, info_cp)
	call2 = "cp %s %s" % (dose, dose_cp) 
	sp.check_call(call1,shell=True)
	sp.check_call(call2,shell=True)
	grm_file = grm_dir+"chr"+str(chrom)
	gcta_call1 = "gcta64 --dosage-mach-gz %s %s --make-grm-gz --out %s" % (dose_cp, info_cp, grm_file)
	cout("Calculating GRM...\n")          
	fout = open("temp_job.sh", 'w')
	script = '''
	#!/bin/bash
	#PBS -N chr%s
	#PBS -S /bin/bash
	#PBS -l mem=10gb
	#PBS -o %sjobs/chr%s.out
	#PBS -e %sjobs/chr%s.err

	module load gcta/1.24.4

	%s 
	rm %s %s   
	''' % (str(chrom),grm_dir,str(chrom),grm_dir,str(chrom),gcta_call1,dose_cp,info_cp)
    
	fout.write(script)
	fout.close()
	cout("Running job....\n")
	job = ["qsub","temp_job.sh"]
	sp.check_call(job)
	os.remove("temp_job.sh")

def make_mgrm():
	global grm_dir
	filename = grm_dir+"grms.txt"
	fout = open(filename,'w')
	for c in range(1,23):
		fout.write(grm_dir+"chr"+str(c)+"\n")
	fout.close()
	gcta_call = "gcta64 --mgrm-gz %s --make-grm-gz --out %s" % (filename,grm_dir+"sc")
	sp.check_call(gcta_call,shell=True) 

def reml():
	global grm_dir, pheno_dir, reml_dir
	str1 = "gcta64 --grm-gz %s --pheno %s" % (grm_dir+"sc", pheno_dir+"sc.pheno")
	str2 = " --qcovar %s --reml --prevalence 0.20 --out %s" % (pheno_dir+"sc.qcovar",reml_dir+"sc") 
	gcta_call = str1 + str2 
	sp.check_call(gcta_call,shell=True)   

def reml_pred():
	global grm_dir, pheno_dir, reml_dir
	str1 = "gcta64 --grm-gz %s --pheno %s" % (grm_dir+"sc", pheno_dir+"sc.pheno")
	str2 = " --qcovar %s --reml --reml-pred-rand --prevalence 0.20 --out %s" % (pheno_dir+"sc.qcovar",reml_dir+"sc") 
	gcta_call = str1 + str2 
	sp.check_call(gcta_call,shell=True)   	
	#return reml_dir+"sc.indi.blp"

def run_blup(chrom):
	'''
	gcta64  --bfile   test   --blup-snp test.indi.blp  --out test
	'''
	global mach_dir, reml_dir, blup_dir 
	dname, iname = "chr"+str(chrom)+".mldose.gz" , "chr"+str(chrom)+".mlinfo.gz" 
	dose, info  = mach_dir+dname, mach_dir+iname 
	dose_cp, info_cp = blup_dir+dname, blup_dir+iname
	cout("Copying dose and info files for chromosome: %s\n" % str(chrom))
	call1 = "cp %s %s" % (info, info_cp)
	call2 = "cp %s %s" % (dose, dose_cp) 
	sp.check_call(call1,shell=True)
	sp.check_call(call2,shell=True)
	blp_file = reml_dir+"sc.indi.blp"
	blup_file = blup_dir+"chr"+str(chrom)
	gcta_call = "gcta64 --dosage-mach-gz %s %s --blup-snp %s --out %s" % (dose_cp, info_cp, blp_file,blup_file)
	cout("Calculating Blup SNP Effects...\n")          
	fout = open("temp_job.sh", 'w')
	script = '''
	#!/bin/bash
	#PBS -N chr%s
	#PBS -S /bin/bash
	#PBS -l mem=10gb
	#PBS -o %sjobs/chr%s.out
	#PBS -e %sjobs/chr%s.err

	module load gcta/1.24.4

	%s 
	rm %s %s   
	''' % (str(chrom),blup_dir,str(chrom),blup_dir,str(chrom),gcta_call,dose_cp,info_cp)
    
	fout.write(script)
	fout.close()
	cout("Running job....\n")
	job = ["qsub","temp_job.sh"]
	sp.check_call(job)
	os.remove("temp_job.sh")

def cat_blup():
	global blup_dir
	fout = gzip.open(blup_dir+"sc.snp.blp.gz",'wb')
	cout("Concatenating BLP SNP files...\n")
	for c in range(1,23):
		cout("Chromosome: %d\n" % c)
		fin = open(blup_dir+"chr"+str(c)+".snp.blp",'r')
		for line in fin:
			l = line.strip()
			fout.write(l+"\n")
	fout.close()

def poly_chrom(chrom):
	'''
	Generate polyscore file for a specificed chromosome 
	'''
	global mach_dir, blup_dir, poly_dir 
	dname, iname = "chr"+str(chrom)+".mldose.gz" , "chr"+str(chrom)+".mlinfo.gz" 
	dose, info  = mach_dir+dname, mach_dir+iname 
	dose_cp, info_cp = poly_dir+dname, poly_dir+iname
	cout("Copying dose and info files for chromosome: %s\n" % str(chrom))
	call1 = "cp %s %s" % (info, info_cp)
	call2 = "cp %s %s" % (dose, dose_cp) 
	sp.check_call(call1,shell=True)
	sp.check_call(call2,shell=True)
	blup_file = blup_dir+"chr"+str(chrom)+".snp.blp"
	out_file = poly_dir+"chr"+str(chrom)+".poly.gz"
	call = "python  %s %s %s %s" % (main_dir+"JTcalcpoly.py",dose_cp,blup_file,out_file)
	cout("Calculating Polyscores...\n")          
	fout = open("temp_job.sh", 'w')
	script = '''
	#!/bin/bash
	#PBS -N chr%s
	#PBS -S /bin/bash
	#PBS -l mem=10gb
	#PBS -o %sjobs/chr%s.out
	#PBS -e %sjobs/chr%s.err

	module load python/2.7.9

	%s 
	rm %s %s   
	''' % (str(chrom),poly_dir,str(chrom),poly_dir,str(chrom),call,dose_cp,info_cp)
    
	fout.write(script)
	fout.close()
	cout("Running job....\n")
	job = ["qsub","temp_job.sh"]
	sp.check_call(job)
	os.remove("temp_job.sh")

def sum_poly():
	global poly_dir
	fout = gzip.open(poly_dir+"sc.poly.gz",'wb')
	cout("Integrating PolyScore Files...\n")
	poly_dic = {}
	call1 = "gunzip %schr1.poly.gz" % poly_dir
	call2 = "gzip %schr1.poly" % poly_dir 
	sp.check_call(call1,shell=True)
	ind_list = list_from_file(poly_dir+"chr1.poly",1)
	sp.check_call(call2,shell=True)
	for ind in ind_list:
		poly_dic[ind] = []  
	for c in range(1,23):
		cout("Chromosome: %d\n" % c)
		fin = gzip.open(poly_dir+"chr"+str(c)+".poly.gz",'rb')
		for line in fin:
			l = line.strip().split()
			iid, score = l[0],l[2]
			poly_dic[iid].append(float(score))
		fin.close()
	for ind in ind_list:
		fout.write(ind+"\t"+ind+"\t"+str(sum(poly_dic[ind]))+"\n")
	#fout.write(l+"\n")
	fout.close()

def make_df_file():
	global pheno_dir, poly_dir
	ind_list = list_from_file(pheno_dir+"sc.pheno",1)
	ind_list.pop(0)
	pdic, qdic, sdic = {}, {}, {}
	fin1 = open(pheno_dir+"sc.pheno",'r')
	for line in fin1:
		l = line.strip().split()
		iid, cs = l[0],l[2]
		pdic[iid] = [cs] 
	fin1.close()
	fin2 = open(pheno_dir+"sc.qcovar",'r')
	for line in fin2:
		l = line.strip().split()
		iid, bmi, age = l[0],l[2],l[3]
		qdic[iid] = [bmi,age]
	fin2.close()
	fin3 = gzip.open(poly_dir+"sc.poly.gz",'rb')
	for line in fin3:
		l = line.strip().split()
		iid, score = l[0],l[2]
		sdic[iid] = [score]
	fin3.close()
	fout = open(pheno_dir+"sc.traits.txt",'w')
	headlist = ["FID","IID","CaseStatus","BMI","AGE","PolyScore"]
	fout.write("\t".join(headlist)+"\n")
	for ind in ind_list:
		write_list = [ind,ind] + pdic[ind] + qdic[ind] + sdic[ind]
		fout.write("\t".join(write_list)+"\n")
	fout.close()



def main():
	#for c in range(1,23):
		#make_grm(c)
	#make_mgrm()
	#reml()
	#reml_pred()
	#for c in range(1,23):
		#run_blup(c)
	#cat_blup()
	#for c in range(1,23):
		#poly_chrom(c)
	#sum_poly()
	make_df_file()


    

if (__name__ == "__main__"): 

	#make_grm(22)
	main() 