#!/public/wangm/miniconda3/bin/python
# -*- coding: utf-8 -*-

import os,time,json
from Bio import SeqIO
from Bio.SeqUtils import GC
from script.config import get_param
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline

param = get_param()
workdir = param[0]
defensefinder = param[3]
blastp = param[4]
blastn = param[5]

tmp_dir = os.path.join(workdir,'tmp')
gb_dir = os.path.join(tmp_dir,'gbk')

VF_Database = os.path.join(workdir,'data','virulence')
IS_Database = os.path.join(workdir,'data','transposase')
arg_Database = os.path.join(workdir,'data','resfinder')
metal_Database = os.path.join(workdir,'data','metal')
pop_Database = os.path.join(workdir,'data','degradation')
sym_Database = os.path.join(workdir,'data','symbiosis')

##### Test
##### Not used
def getmmseq(runID):

	mmseq_dir = os.path.join(tmp_dir,runID,'mmseq')
	if not os.path.exists(mmseq_dir):
		os.makedirs(mmseq_dir)
	
	if not os.path.exists(os.path.join(tmp_dir,runID,runID+'.locus_tag.faa')):
		infaa = os.path.join(gb_dir,runID+'.faa')
	else:
		infaa = os.path.join(tmp_dir,runID,runID+'.locus_tag.faa')

	argDB = os.path.join('mmseqDB','resfinder','resfinder')
	vfDB = os.path.join('mmseqDB','vfdb','vfdb')
	isDB = os.path.join('mmseqDB','isfinder','isfinder')

	tmpd = os.path.join(mmseq_dir,'tmp')
	queryDB = os.path.join(mmseq_dir,runID+'_queryDB')
	argres = os.path.join(mmseq_dir,runID+'_tarDB','resfDB')
	vfres = os.path.join(mmseq_dir,runID+'_tarDB','vfresDB')
	isres = os.path.join(mmseq_dir,runID+'_tarDB','isresDB')
	argout = os.path.join(tmp_dir,runID,'arg.m8')
	vfout = os.path.join(tmp_dir,runID,'vf.m8')
	isout = os.path.join(tmp_dir,runID,'is.m8')

	if not os.path.exists(os.path.join(mmseq_dir,runID+'_tarDB')):
		os.makedirs(os.path.join(mmseq_dir,runID+'_tarDB'))	

	queryDBc = ['mmseqs','createdb',infaa,queryDB, '> /dev/null']
	argc = ['mmseqs', 'search', queryDB, argDB, argres, tmpd, '> /dev/null']
	argconv = ['mmseqs', 'convertalis', '--format-mode 2', queryDB, argDB, argres, argout, '> /dev/null']

	vfc = ['mmseqs', 'search', queryDB, vfDB, vfres, tmpd, '> /dev/null']
	vfconv = ['mmseqs', 'convertalis', '--format-mode 2', queryDB, vfDB, vfres, vfout, '> /dev/null']

	isc = ['mmseqs', 'search', queryDB, isDB, isres, tmpd, '> /dev/null']
	isconv = ['mmseqs', 'convertalis', '--format-mode 2', queryDB, isDB, isres, isout, '> /dev/null']

	os.system(' '.join(queryDBc))
	os.system(' '.join(argc))
	os.system(' '.join(vfc))
	os.system(' '.join(isc))
	os.system(' '.join(argconv))
	os.system(' '.join(vfconv))
	os.system(' '.join(isconv))

	argdict = getres(argout)
	vfdict = getres(vfout)
	isdict = getres(isout)
	dfdict = getdf(runID)

	return argdict,vfdict,isdict,dfdict

def getres(mmout): ##not used

	mmdict = {}
	with open(mmout) as mmseqr:
		for line in mmseqr.readlines():
			lines = line.strip().split('\t')
			hvalue = int(lines[3])/int(lines[13])*float(lines[2])
			if lines[0] not in mmdict:
				if hvalue >= 0.64:
					mmdict[lines[0]] = lines[1].split('|')[1]
	return mmdict

#### Not used
#### Test

def getdf(runID):

	if not os.path.exists(os.path.join(tmp_dir,runID,runID+'.locus_tag.faa')):
		infaa = os.path.join(gb_dir,runID+'.faa')
	else:
		infaa = os.path.join(tmp_dir,runID,runID+'.locus_tag.faa')

	dfout = os.path.join(tmp_dir,runID,'defense_'+runID)
	defcmd = [defensefinder, 'run', '-w 8 --models-dir ./data/macsydata/','-o', dfout, infaa, '> /dev/null']
	os.system(' '.join(defcmd))

	dfdict = {}
	with open(os.path.join(dfout,'defense_finder_genes.tsv')) as dfres:
		for line in dfres.readlines():
			lines = line.strip().split('\t')
			if lines[0] != 'replicon':
				dfdict[lines[1]] = lines[2].replace('__',',')

	return dfdict

def isblast(faa_file,IS_out):
	blastp_cline = NcbiblastpCommandline(cmd=blastp, query=faa_file, db=IS_Database, \
                       evalue=0.0001, num_threads=20, max_hsps=1, num_descriptions=1, \
                       num_alignments=1, outfmt="6 std slen stitle", out=IS_out)
	blastp_cline()

def vfblast(faa_file,VF_out):
	blastp_cline = NcbiblastpCommandline(cmd=blastp, query=faa_file, db=VF_Database, \
                       evalue=0.0001, num_threads=20, max_hsps=1, num_descriptions=1, \
                       num_alignments=1, outfmt="6 std slen stitle", out=VF_out)
	blastp_cline()

def argblast(fa_file,arg_out):
	blastp_cline = NcbiblastnCommandline(cmd=blastn, query=fa_file, db=arg_Database, \
                       evalue=0.0001, num_threads=20, max_hsps=1, num_descriptions=1, \
                       num_alignments=1, outfmt="6 std slen stitle", out=arg_out)
	blastp_cline()

def metalblast(faa_file,metal_out):
	blastp_cline = NcbiblastpCommandline(cmd=blastp, query=faa_file, db=metal_Database, \
                       evalue=0.0001, num_threads=20, max_hsps=1, num_descriptions=1, \
                       num_alignments=1, outfmt="6 std slen stitle", out=metal_out)
	blastp_cline()

def popblast(faa_file,pop_out):
	blastp_cline = NcbiblastpCommandline(cmd=blastp, query=faa_file, db=pop_Database, \
                       evalue=0.0001, num_threads=20, max_hsps=1, num_descriptions=1, \
                       num_alignments=1, outfmt="6 std slen stitle", out=pop_out)
	blastp_cline()

def symblast(faa_file,sym_out):
	blastp_cline = NcbiblastpCommandline(cmd=blastp, query=faa_file, db=sym_Database, \
                       evalue=0.0001, num_threads=20, max_hsps=1, num_descriptions=1, \
                       num_alignments=1, outfmt="6 std slen stitle", out=sym_out)
	blastp_cline()

def havalue(value,out):

	blast_filter = {}
	for line in open(out,'r').readlines():
		lines = line.strip().split('\t')
		havalue = (int(lines[3])/int(lines[12]))*float(lines[2])/100
		if havalue >= float(value):
			blast_filter[lines[0]]=lines[1].split('|')[1]
	return blast_filter

def getblast(runID):
	
	arg_out = os.path.join(tmp_dir,runID,'arg.m8')
	vf_out = os.path.join(tmp_dir,runID,'vf.m8')
	is_out = os.path.join(tmp_dir,runID,'is.m8')
	pop_out = os.path.join(tmp_dir,runID,'pop.m8')
	metal_out = os.path.join(tmp_dir,runID,'metal.m8')
	sym_out = os.path.join(tmp_dir,runID,'sym.m8')


	if not os.path.exists(os.path.join(tmp_dir,runID,runID+'.locus_tag.faa')):
		infaa = os.path.join(gb_dir,runID+'.faa')
		infa = os.path.join(gb_dir,runID+'.ffn')
	else:
		infaa = os.path.join(tmp_dir,runID,runID+'.locus_tag.faa')
		infa = os.path.join(tmp_dir,runID,runID+'.locus_tag.spaceHeader.ffn')

	isblast(infaa,is_out)
	vfblast(infaa,vf_out)
	argblast(infa,arg_out)

	metalblast(infaa,metal_out)
	popblast(infaa,pop_out)
	symblast(infa,sym_out)

	isdict = havalue('0.64',is_out)
	vfdict = havalue('0.64',vf_out)
	argdict = havalue('0.81',arg_out)
	metaldict = havalue('0.64',metal_out)
	popdict = havalue('0.64',pop_out)
	symdict = havalue('0.64',sym_out)

	dfdict = getdf(runID)

	return argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict




