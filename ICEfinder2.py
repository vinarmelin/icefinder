#!/usr/bin/env python3

# -*- coding: utf-8 -*-
#
# ICEfinder: Detecting Integrative and Conjugative Elements in Bacteria.
# Meng Wang on Sep-4-2023
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University
# Version 2.0 - Sep 4, 2023
########################################################################

import os
import argparse
import logging
from script.checkin import get_fagb
from script.single import _single
from script.metaICE import _meta
from script.config import get_param

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

workdir,kraken,krakenDB,defensefinder,blastp,blastn,seqkit,prodigal,prokka,macsyfinder,hmmsearch = get_param()
logging.info(f"""
###################
Starting ICEfinder2
###################

config.ini data -->
	Workdir: {workdir}
	kraken2: {kraken}
	kraken2 DB: {krakenDB}
	defensefinder: {defensefinder}
	blastp: {blastp}
	blastn: {blastn}
	seqkit: {seqkit}
	prodigal: {prodigal}
	prokka: {prokka}
	macsyfinder: {macsyfinder}
	hmmsearch: {hmmsearch}
""")

tmp_dir = os.path.join(workdir,'tmp')
fa_dir = os.path.join(tmp_dir,'fasta') 
gb_dir = os.path.join(tmp_dir,'gbk')

def add_arguments_to_parser(parser):

	parser.add_argument('-v', '--version', action='version', version='2.0',
                        help="Show ICEfinder version")
	parser.add_argument('-i', '--input', type=str, required=True,
                        help='FASTA/Genbank format file, Genbank format file accepted only for single genome.')
	parser.add_argument('-t', '--type', type=str, required=True,
                        help='Genome Type: Single/Metagenome')

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='ICEfinder', usage='python ICEfinder.py -i fasta_file/genbank_file -t Single/Metagenome',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	add_arguments_to_parser(parser)
	args = parser.parse_args()
	intype = args.type
	input_file = args.input

	if not os.path.exists(tmp_dir):
		os.mkdir(tmp_dir)
	if not os.path.exists(gb_dir):
		os.mkdir(gb_dir)
	if not os.path.exists(fa_dir):
		os.mkdir(fa_dir)

	runID = file_name_without_extension = os.path.splitext(os.path.basename(input_file))[0]

	logging.info(f"""
Input data -->
	File: {input_file}
	intype: {intype}
	runID: {runID}
	Tmp folder: {tmp_dir}
	Genebank folder: {gb_dir}
	Fasta folder: {fa_dir}
""")

	infile,filetype = get_fagb(runID,input_file,intype)

	if intype == 'Single':
		logging.info("Executing Single mode")
		_single(runID,infile,filetype)
	else:
		logging.info("Executing Meta mode")
		_meta(runID,infile)

	print(runID+' done!!')
