#!/public/wangm/miniconda3/bin/python
# -*- coding: utf-8 -*-

import os,sys,shutil
from Bio import SeqIO
from script.config import get_param

param = get_param()
workdir = param[0]
tmp_dir = os.path.join(workdir,'tmp') 
in_dir = os.path.join(tmp_dir,'fasta')
gb_dir = os.path.join(tmp_dir,'gbk')

def is_fagb(filename):

	filetype = ''
	with open(filename, "r") as handle1:
		fasta = SeqIO.parse(handle1, "fasta")
		if any(fasta):
			filetype = 'fa'
	with open(filename, "r") as handle2:
		gbk = SeqIO.parse(handle2, "gb")
		if any(gbk):
			filetype = 'gb'
	return filetype

def remove_folders_with_runID(root_dir, runID):
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
        for dirname in dirnames:
            if runID in dirname:
                dir_to_remove = os.path.join(dirpath, dirname)
#                print(f"Removing folder: {dir_to_remove}")
                shutil.rmtree(dir_to_remove)
        for filename in filenames:
            if runID in filename:
                file_to_remove = os.path.join(dirpath, filename)
#                print(f"Removing file: {file_to_remove}")
                os.remove(file_to_remove)

def get_fagb(runID,input_file,intype):

#	for filename in os.listdir(gb_dir):
#		file_path = os.path.join(gb_dir, filename)
#		if runID in filename:
#			os.remove(file_path)

#	folder_path = os.path.join(tmp_dir, runID)
#	if os.path.exists(folder_path) and os.path.isdir(folder_path):
#		shutil.rmtree(folder_path)
	remove_folders_with_runID(gb_dir,runID)
	remove_folders_with_runID(tmp_dir,runID)

	infile = os.path.join(in_dir,runID)
	shutil.copy(input_file, infile)
	try:
		filetype = is_fagb(infile)
	except:
		print('ERROR: The input file is not a standard FASTA/GenBank format! Please check !')
		sys.exit()

	if not filetype:
		print('ERROR: The input file is not a standard FASTA/GenBank format! Please check !')
		sys.exit()

	else:
		if filetype == 'fa':
			newfile = os.path.join(in_dir,runID+'.fa')
			seq_record = SeqIO.parse(infile, 'fasta')
			if len(list(seq_record)) == 1:
				for seq_records in SeqIO.parse(infile, 'fasta'):
					seq_records.id = runID
					SeqIO.write(seq_records, newfile, 'fasta')
			elif intype == 'Metagenome':
				shutil.copy(infile, newfile)
				filetype = 'multifa'
			else:
				print('ERROR: Input file accepted for one sequence only.')
				sys.exit()
		else:
			if not os.path.exists(gb_dir):
				os.makedirs(gb_dir)
				os.chmod(gb_dir, 0o777)
			gbfile = os.path.join(gb_dir,runID+'.gbk')
			newfile = os.path.join(in_dir,runID+'.fa')
			seq_record = SeqIO.parse(infile, 'gb')
			if len(list(seq_record)) == 1:
				i = 0
				for seq_records in SeqIO.parse(infile, 'gb'):
					for seq_feature in seq_records.features:
						if seq_feature.type=="CDS":
							if 'locus_tag' in seq_feature.qualifiers:
								continue
							else:
								i += 1
								if i > 10:
									print('ERROR: Too many CDS do not have locus_tag in GenBank input file! Please check or try FASTA format input!')
									sys.exit()
					if len(set(list(str(seq_records.seq)))) == 1:
						print('ERROR: The uploaded file is not a standard GenBank format! Please check or try a FASTA format input!')
						sys.exit()

					else:			
						seq_records.id = runID
						SeqIO.write(seq_records, gbfile, 'gb')
						SeqIO.write(seq_records, newfile, 'fasta')
			else:
				print('ERROR: Input file accepted for one sequence only.')
				sys.exit()

	return newfile,filetype
