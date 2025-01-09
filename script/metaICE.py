#!/public/wangm/miniconda3/bin/python
# -*- coding: utf-8 -*-

import os,time
import random,json
import string,shutil
from Bio import SeqIO
from ete3 import NCBITaxa
from Bio.SeqUtils import GC
from functools import cmp_to_key
from script.function import getblast
from script.config import get_param

param = get_param()
workdir = param[0]
soft_dir = param[1]
kraken = param[1]
krakenDB = param[2]
blastn = param[5]
seqkit = param[6]
prodigal = param[7]
prokka = param[8]
macsyfinder = param[9]
hmmsearch = param[10]

tmp_dir = os.path.join(workdir,'tmp') 
in_dir = os.path.join(tmp_dir,'fasta')
gb_dir = os.path.join(tmp_dir,'gbk')

def rename(runID, infile):
	run_dir = os.path.join(tmp_dir,runID)
	if not os.path.exists(run_dir):
		os.makedirs(run_dir)
	if not os.path.exists(gb_dir):
		os.mkdir(gb_dir)

	filename = os.path.basename(infile)
	resultf = filename.rsplit('.', 1)[0]
	infile1 = os.path.join(os.path.dirname(infile), resultf)

	i = 1
	id_dict = {}
	newIDfa = os.path.join(run_dir,runID+'_newID.fa')
	newfa = open(newIDfa,'w')

	for seq_record in SeqIO.parse(infile, "fasta"):
		if len(seq_record.id) > 15:
			realID = seq_record.id[:15]+'...'
		else:
			realID = seq_record.id
		contigID = 'contig_' + str(i)
		seq_record.id = contigID
		seqfa = str(seq_record.seq)
		newfa.write(">%s\n%s\n" % (contigID,seqfa))
		id_dict[contigID] = realID
		i += 1
	return id_dict

def get_time():

	return time.asctime( time.localtime(time.time()) )

def Taxonomy(runID):

	newIDfa = os.path.join(tmp_dir,runID,runID+'_newID.fa')
	report = os.path.join(tmp_dir,runID,runID+'_kraken.report')
	output = os.path.join(tmp_dir,runID,runID+'_kraken.output')
	drawout = os.path.join(tmp_dir,runID,'kraken.html')

	annocmd = "{} {} {} {} {} {} {} {} {}".format(kraken,"--db",krakenDB,"--report",\
		report,"--output",output,newIDfa,'>/dev/null')
	os.system(annocmd)
#	drawcmd = ' '.join(['/opt/R/3.6.3/bin/Rscript', './script/sankey.R', report, drawout, '>/dev/null'])
#	os.system(drawcmd)

	spdict = {}
	with open(output,'r') as taxainfo:
		for line in taxainfo.readlines():
			lines = line.strip().split('\t')
			ID,taxid = lines[1],lines[2]
			if taxid == '0':
				spdict[ID] = '-'
			else:
				spname,strainame = get_ranks(taxid)
				spdict[ID] = spname

	return drawout,spdict,report

def get_ranks(taxid):
	ncbi = NCBITaxa()
	lineage = ncbi.get_lineage(taxid)   
	names = ncbi.get_taxid_translator(lineage)
	lineage2ranks = ncbi.get_rank(names)
	ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())

	strainid = ''
	if 'species' in ranks2lineage:
		spid = ranks2lineage['species']
		if 'strain' in ranks2lineage:
			strainid = ranks2lineage['strain']
		else:
			strainid = ''
	elif 'genus' in ranks2lineage:
		spid = ranks2lineage['genus']
		strainid = ''
	elif 'phylum' in ranks2lineage:
		spid = ranks2lineage['phylum']
		strainid = ''

	try:
		spname = list(ncbi.get_taxid_translator([spid]).values())[0]
	except:
		spname = '-'

	strainame = ''
	if strainid:
		strainame = list(ncbi.get_taxid_translator([strainid]).values())[0]

	return spname,strainame

def getbase(runID):

	newIDfa = os.path.join(tmp_dir, runID, runID+'_newID.fa')
	statcmd = ' '.join([seqkit,"stats","-a",newIDfa])
	stats = os.popen(statcmd, "r")
	for line in stats.readlines():
		lines = line.strip().split()
		if lines[0] != 'file':
			lengt = lines[4]
			count = lines[3]
			n50 = lines[12]

	basefile = os.path.join(tmp_dir,runID, runID+'_info.json')
	basedict = {'JobID':runID,
		    'Submission date':get_time(),
		    'Total length': lengt+' bp',
   		    'Contig number': count,
		    'Sequence N50': n50+' bp'
		    }
	with open(basefile,'w') as gein2:
		json.dump(basedict, gein2, indent=4)

	return basefile

def gc(fasta_file,start,end):

	record = SeqIO.read(fasta_file, "fasta")
	sequence = record.seq[start-1:end]
	gcs = str("%.2f"%GC(sequence))

	return gcs

def calculate_gc(fasta_file, start, end, window_size, step_size):

	record = SeqIO.read(fasta_file, "fasta")
	if start == 0:
		start = 1
	sequence = record.seq[start-1:end]

	windows = []
	gc_contents = []
	pos = []
	j = start/1000 + 0.025
	for i in range(0, len(sequence) - window_size + 1, step_size):
		window = sequence[i:i+window_size]
		gc_content = GC(window)
		gc_contents.append(gc_content)
		pos.append(round(j, 4))
		j += 0.05

	gcdict = {
		        'xData':pos,
		        'datasets':[{
		            'name':'',
		            'data':gc_contents,
		            'unit':'%',
		            'type':'line',
		            "valueDecimals": 1
		        }]
	    }

	return gcdict

def preanno(runID):

	newIDfa = os.path.join(tmp_dir, runID, runID+'_newID.fa')
	anno_fa = os.path.join(tmp_dir, runID, runID + '.faa')
	anno_gff = os.path.join(tmp_dir, runID, runID + '.gff')
	anno_cmd = prodigal +' -c -m -q -p meta -f gff -i ' + newIDfa + ' -a ' + anno_fa + ' -o ' + anno_gff
	os.system(anno_cmd)

def scanf(hmmlist):

	ICEcount = []
	for line in hmmlist:
		if 'MOB' in line:
			ICEcount.append('MOB')
		elif 't4cp' in line or 'tcpA' in line:
			ICEcount.append('t4cp')
		elif 'FA' in line:
			ICEcount.append('T4SS')
		elif line in ['Phage_integrase','UPF0236','Recombinase','rve','TIGR02224','TIGR02249','TIGR02225','PB001819']:
			ICEcount.append('Int')
		else:
			ICEcount.append('T4SS')
	if ICEcount.count('MOB') and ICEcount.count('t4cp') and ICEcount.count('Int') and ICEcount.count('T4SS') >= 5:
		return True
	else:
		return False

def prescan(runID):
	preanno(runID)
	anno_fa = os.path.join(tmp_dir, runID, runID + '.faa')
	scanfile = os.path.join(tmp_dir, runID, runID + '_prescan')
	scancmd = [ './tool/hmmscan2', '--tblout', scanfile, './data/ICEscan.hmm', anno_fa, '> /dev/null']
	os.system(' '.join(scancmd))

	icedict = {}
	chosen = []
	with open(scanfile, 'r') as outfile:
		for line in outfile.readlines():
			if not line.startswith('#'):
				lines = line.strip().split()
				if lines[2] in icedict:
					continue
				id_parts = lines[2].split('_')
				key = '_'.join(id_parts[0:2])
				if float(lines[4]) < 0.00001:
					if key in icedict:
						icedict[key].append(lines[0])
					else:
						icedict[key] = [lines[0]]
	for k,v in icedict.items():
		if scanf(v):
			chosen.append(k)

	return chosen

def prokkanno(runID,infile):

	cmd = [prokka,infile,"--force","--fast --quiet --cdsrnaolap --cpus 8 ","--outdir",gb_dir,"--prefix",runID]
	os.system(' '.join(cmd))

def ICEscan(runID):

	anno_fa = os.path.join(gb_dir, runID + '.faa')
	ICE_res = os.path.join(tmp_dir, runID, runID + '_ICE')
	ICE_cmd = macsyfinder + ' --db-type ordered_replicon --hmmer '+ hmmsearch + ' --models-dir ./data/macsydata/ --models ICEscan all --replicon-topology linear --coverage-profile 0.3 --sequence-db '+ anno_fa+ ' -o ' + ICE_res + ' > /dev/null'
	os.system(ICE_cmd)

def getgff(runID):

	gffile = os.path.join(gb_dir, runID + '.gff')
	with open(gffile,'r') as gffin:
		trnadict = {}
		posdict = {}
		for line in gffin.readlines():
			if 'ID=' in line:
				lines = line.strip().split('\t')
				ids = lines[8].split(';')[0].split('=')[1]
				header = ids.split('_')[0]
				product = lines[8].split('product=')[1]
				pos = [lines[3],lines[4],lines[6],product]
				if lines[2] == 'tRNA' or lines[2] == 'tmRNA':
					trnadict[ids] = product
				posdict[ids] = pos
				totalnum = getnum(ids)

	return trnadict,posdict,header,totalnum

def getnum(ID):

	return int(ID.split('_')[1].lstrip("0"))

def find_max_distance(numbers):
	max_distance = -1
	max_distance_index = -1

	for i in range(len(numbers) - 1):
		distance = abs(numbers[i] - numbers[i+1])
		if distance > max_distance:
			max_distance = distance
			max_distance_index = i

	if max_distance_index == -1:
		return None

	return numbers[max_distance_index], numbers[max_distance_index + 1]

def pos_tag(pos,posdict,ICE,final,totalnum,dirtag):

	tICE = ICE
	tfinal = final
	for k,v in posdict.items():
		vstart,vend = int(v[0]),int(v[1])
		if int(pos) <= vend:
			if dirtag == 's':
				tICE = getnum(k)
				tfinal = max(1, tICE - 5)
			else:
				if vstart > int(pos):
					tICE = getnum(k) - 1 
				else:
					tICE = getnum(k)
				tfinal = min(totalnum, tICE + 5)
			break
	return tICE, tfinal

def merge_tRNA(runID,ICEdict,DRlist):

	trnadict,posdict,header,totalnum = getgff(runID)
	fICE = getnum(next(iter(ICEdict)))
	eICE = getnum(list(ICEdict.keys())[-1])

	nfICEnum = max(1, fICE - 5)
	neICEnum = min(totalnum, eICE + 5)

	ICEtagnum = [nfICEnum,neICEnum]
	trnalist = []
	for key,value in trnadict.items():
		if nfICEnum <= getnum(key) <= neICEnum:
			ICEtagnum.append(getnum(key))
			trnalist.append(value)

	ICEtagnum.sort()
	finalstart,finalend = find_max_distance(ICEtagnum)

	myDR1 = posdict[zill(header,fICE)][0]
	myDR2 = ''
	myDR3 = ''
	myDR4 = posdict[zill(header,eICE)][1]

	if trnalist:
		if finalstart == nfICEnum:
			eICE = finalend
			finalend = min(totalnum, finalend + 5)
			myDR4 = posdict[zill(header,eICE)][1]					
			for line in DRlist:
				DRs = line.split('|')
				if int(DRs[3]) - int(DRs[0]) > 500000:
					continue
				if int(DRs[3]) - int(DRs[0]) < 5000:
					continue					
				if int(posdict[zill(header,eICE)][0]) < int(DRs[3]) < int(posdict[zill(header,eICE)][1]):
					checktrna = 0
					for key,value in trnadict.items():
						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break

					fICE,finalstart = pos_tag(DRs[0],posdict,fICE,finalstart,totalnum,'s')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]
					break

		elif finalend == neICEnum:
			fICE = finalstart
			finalstart =  max(1, finalstart - 5)
			myDR1 = posdict[zill(header,fICE)][0]
			for line in DRlist:
				DRs = line.split('|')
				if int(DRs[3]) - int(DRs[0]) > 500000:
					continue	
				if int(DRs[3]) - int(DRs[0]) < 5000:
					continue									
				if int(posdict[zill(header,fICE)][0]) < int(DRs[0]) < int(posdict[zill(header,fICE)][1]):
					checktrna = 0
					for key,value in trnadict.items():
						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break
					eICE,finalend = pos_tag(DRs[3],posdict,eICE,finalend,totalnum,'e')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]									
					break

	return myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist

def get_DR(runID,infile):

	DRindex = os.path.join(tmp_dir, runID, runID+'_DR')
	DRout = os.path.join(tmp_dir, runID, runID+'_DRout')
	maktree_cmd = './tool/mkvtree -db ' + infile +' -indexname '+ DRindex +' -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1'
	vmatch_cmd = './tool/vmatch -l 15 '+ DRindex + ' > '+ DRout
	os.system(maktree_cmd)
	os.system(vmatch_cmd)

	DRlist = []
	with open(DRout,'r') as DRin:
		for line in DRin.readlines():
			lines = line.strip().split()
			if not line.startswith('#'):
				DR = [lines[2],str(int(lines[2])+int(lines[0])),lines[6],str(int(lines[6])+int(lines[4]))]
				DRlist.append('|'.join(DR))
	return DRlist

def get_ICE(runID,infile):

	ICE_dir = os.path.join(tmp_dir, runID, runID + '_ICE')
	ICE_res = os.path.join(ICE_dir,'all_systems.tsv')

	if os.path.exists(ICE_dir):
		os.system('rm -r ' + ICE_dir)
	ICEscan(runID)

	with open(ICE_res,'r') as ICEin:
		ICEdict = {}
		infodict = {}
		for line in ICEin.readlines():

			if 'Chromosome' in line:
				lines = line.strip().split('\t')
				if 'UserReplicon_IME' not in lines[5]:
					gbname = lines[1]
					tags = get_feat(lines[2])
					mpf = lines[4].split('/')[-1].split('_')[1]

					if 'Relaxase@' in tags:
						mob = tags.split('@')[1]
					else:
						mob = ''
					ICEtag = 'ICE'+lines[5].split('_')[-1]

					ICEdict.setdefault(ICEtag,{})[gbname]=tags

					if ICEtag not in infodict:
						infodict[ICEtag] = {'mob': [], 'mpf': []}
					if mob not in infodict[ICEtag]['mob']:
						if mob:
							infodict[ICEtag]['mob'].append(mob)
					if mpf not in infodict[ICEtag]['mpf']:
						infodict[ICEtag]['mpf'].append(mpf)	

	dictICE = {}
	posdict = {}
	trnalist = []
	header = ''
	DRlist = get_DR(runID,infile)
	for key,value in ICEdict.items():
		myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist = merge_tRNA(runID,value,DRlist)
		dictICE[key] = [myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend]

	return dictICE,ICEdict,posdict,header,trnalist,infodict

def args(runID):

	return getblast(runID)

def zill(header,num):

	return header+ '_' +str(num).zfill(5)

def get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product):

	feature = [feature]
	product = [product]

	if gene in argdict:
		feature.append('AR')
		product.append(argdict[gene])
	if gene in vfdict:
		feature.append('VF')
		product.append(vfdict[gene])
	if gene in isdict:
		feature.append('IS')
		product.append(isdict[gene])
	if gene in dfdict:
		feature.append('Defense')
		product.append(dfdict[gene])

	if gene in metaldict:
		feature.append('Metal')
		product.append(metaldict[gene])
	if gene in popdict:
		feature.append('Degradation')
		product.append(popdict[gene])
	if gene in symdict:
		feature.append('Symbiosis')
		product.append(symdict[gene])

	feature = '; '.join(list(filter(None, feature)))
	product = '; '.join(list(filter(None, product)))

	return feature,product

def oritseq(runID, regi, infile, start, end):

	oritseq = '-'
	fafile = os.path.join(tmp_dir,runID,regi+'_fororit.fa')	
	with open(fafile,'w') as orif:
		seq = getfa(infile,start,end)
		orif.write('>fororit\n')
		orif.write(seq)

	oriT_Database = os.path.join(workdir,'data','oriT_db')
	blastn_out = os.path.join(tmp_dir,runID,regi+'_oriTout')
	blast_cmd = [blastn, "-db", oriT_Database, "-query", fafile, "-evalue 0.01 -word_size 11 -outfmt '6 std qlen slen' -num_alignments 1 -out", blastn_out,">/dev/null"]
	os.system(' '.join(blast_cmd))

	with open(blastn_out,'r') as oritout:
		for line in oritout.readlines():
			lines = line.strip().split()
			if lines[0]:
				matchl = int(lines[3])
				slen = int(lines[13])
				ident = float(lines[2])
				hvalue = (matchl/slen)*ident	
				if hvalue > 0.49:
					oritseq = getfa(fafile,str(int(lines[6])-1),lines[7])
					break
	return oritseq

def get_feat(feat):

	featuredict = {
		'Phage_integrase':'Integrase','UPF0236':'Integrase',
        'Recombinase':'Integrase','rve':'Integrase',
        'TIGR02224':'Integrase','TIGR02249':'Integrase',
        'TIGR02225':'Integrase','PB001819':'Integrase'
	}
	if feat in featuredict:
		return 'Integrase@'+feat
	elif 'T4SS_MOB' in feat:
		tag = feat.split('_')[1]
		return 'Relaxase@'+tag
	elif 't4cp' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	elif 'tcpA' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	else:
		return 'T4SS@'+feat.replace('T4SS_','')

def getcolor(feature,product):

	coldict = {'DR':'black','Gene':'#C0C0C0',
		   'Hyp':'#DCDCDC','Integrase':'blue',
		   'Transposase':'yellow','T4SS':'lightpink','T4CP':'orange',
		   'Relaxase':'brown','AR':'red','tRNA':'black',
		   'Flank':'gray','VF':'#ba8448','Defense':'#00B050',
		   'Metal':'#03A89E','Degradation':'#640B0F','Symbiosis':'#FFFFCD'
	}

	namedict = {'Hyp':'Hypothetical protein','Gene':'Other gene',
		    'AR':'Antibiotic resistance gene',
		    'VF':'Virulence factor','Metal':'Metal resistance',
		    'Flank':'Flank region','Defense':'Defense system',
		    'Transposase':'Transposase','Relaxase':'Relaxase',
		    'T4CP':'T4CP','T4SS':'T4SS','Integrase':'Integrase',
		    'Degradation':'Degradation','Symbiosis':'Symbiosis'
	}

	if 'Integrase' in feature:
		feature = 'Integrase'
	elif 'T4SS' in feature:
		feature = 'T4SS'
	elif 'T4CP' in feature:
		feature = 'T4CP'
	elif 'Relaxase' in feature:
		feature = 'Relaxase'
	elif 'IS' in feature:
		feature = 'Transposase'
	elif 'VF' in feature:
		feature = 'VF'
	elif 'AR' in feature:
		feature = 'AR'
	elif 'Defense' in feature:
		feature = 'Defense'
	elif 'Metal' in feature:
		feature = 'Metal'
	elif 'Degradation' in feature:
		feature = 'Degradation'
	elif 'Symbiosis' in feature:
		feature = 'Symbiosis'

	elif feature == 'Flank':
		feature == 'Flank'
	elif feature == '':
		if product == 'hypothetical protein':
			feature = 'Hyp'
		else:
			feature = 'Gene'
	else:
		feature = 'Gene'

	return coldict[feature], namedict[feature]

def gstrand(instra):

	strands = {'+' : 1, '-' : -1}
	return strands[instra]

def getfa(infile,s,e):

	seq_record = SeqIO.read(infile, "fasta")
	sequence = seq_record.seq[int(s):int(e)]
	return str(sequence)

def get_map(sprunID,spdict,id_dict):

	final_dir = os.path.join(workdir,'tmp',sprunID,'result')
	js_dir = os.path.join(final_dir,'js')
	if not os.path.exists(js_dir): 
		os.makedirs(js_dir)
	gcmap = os.path.join(workdir,'script','js','gcmap.js')
	viewfile = os.path.join(workdir,'script','js','view.html')

	fasta_file = os.path.join(tmp_dir, sprunID, sprunID+'.fa')
	dictICE,ICEdict,posdict,header,trnalist,infodict = get_ICE(sprunID,fasta_file)
	argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict = args(sprunID)

	ICEss = {}
	for key,value in dictICE.items():
		genelist = []
		regi = sprunID+'_'+key
		regijs = 'contig_'+sprunID.split("_contig_", 1)[-1] +'_'+key
		genefile = os.path.join(final_dir,regi+'_gene.json')
		infofile = os.path.join(final_dir,regi+'_info.json')
		gcjson = os.path.join(js_dir,regijs+'_gc.js')
		mapfile = os.path.join(js_dir,regijs+'.js')
		htmlfile = os.path.join(final_dir,regi+'.html')
		[myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend] = value

		start = finalstart
		while start < fICE:
			gene = zill(header,start)
			s,e,strand,pro = posdict[gene]
			pos = s+'..'+e+' ['+strand+'], '+str(int(e)-int(s)+1)

			feature = 'Flank'
			product = pro
			feature,product = get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product)
			if 'hypothetical protein;' in product:
				product = product.replace('hypothetical protein;','')

			start += 1
			content = {
					'gene':gene,
					'pos':pos,
					'prod': product,
					'featu': feature
			}
			genelist.append(content)

		mov = fICE
		while mov <= eICE:
			gene = zill(header,mov)
			s,e,strand,pro = posdict[gene]
			pos = s+'..'+e+' ['+strand+'], '+str(int(e)-int(s)+1)

			if gene in ICEdict[key]:
				[feature,pro11] = ICEdict[key][gene].split('@')
			else:
				feature,pro11 = '',''

			if pro11:
				if pro == 'hypothetical protein':
					product = pro11
				else:
					product = pro+', '+ pro11
			else:
				product = pro

			feature,product = get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product)
			mov += 1
			content = {
					'gene':gene,
					'pos':pos,
					'prod': product,
					'featu': feature
			}
			genelist.append(content)				

		while mov <= finalend:
			gene = zill(header,mov)
			s,e,strand,pro = posdict[gene]
			pos = s+'..'+e+' ['+strand+'], '+str(int(e)-int(s)+1)

			feature = 'Flank'
			product = pro
			feature,product = get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product)
			if 'hypothetical protein;' in product:
				product = product.replace('hypothetical protein;','')

			mov += 1
			content = {
					'gene':gene,
					'pos':pos,
					'prod': product,
					'featu': feature
			}
			genelist.append(content)

		with open(genefile,'w') as gene_file:
			json.dump(genelist, gene_file, indent=4)

		contigID = sprunID.split('_', 1)[1]

		sgene = zill(header,fICE)
		egene = zill(header,eICE)
		s1,e1,strand1,pro1 = posdict[sgene]
		s2,e2,strand2,pro2 = posdict[egene]
		if myDR1 == '0':
			myDR1 = '1'

		ICEss[regi] = '|'.join([myDR1,myDR4,str(fICE),str(eICE)])

		host = spdict[contigID]
		gcc = gc(fasta_file,int(myDR1),int(myDR4))
		source = id_dict[contigID]

		if myDR2:
			DR1 = getfa(fasta_file,myDR1,myDR2)
			DR2 = getfa(fasta_file,myDR3,myDR4)
			DRw = 'attL:'+myDR1+'..'+myDR2+'('+DR1+')  '+'attR:'+myDR3+'..'+myDR4+'('+DR2+')'
		else:
			DRw = '-'

		oritseqs = oritseq(sprunID, regi, fasta_file, myDR1,myDR4)
#		oritdesc = "<br>".join([oritseqs[i:i+63] for i in range(0, len(oritseqs), 63)])

		ICEinfo = {
			'Contig source':source,
			'Host Strain':host,
			'GC Content (%)':gcc,
			'Length (bp)':str(int(e2)-int(s1)+1),
			'oriT seq':oritseqs,
			'DRs':DRw,
			'Relaxase Type': ','.join(infodict[key]['mob']),
			'Mating pair formation systems':','.join(infodict[key]['mpf']),
			'Close to tRNA':','.join(trnalist)
		}
		with open(infofile,'w') as info_file:
			json.dump(ICEinfo, info_file, indent=4)

		i = 1
		mapzlist = []
		mapflist = []
		for gene in genelist:
			color, name = getcolor(gene['featu'],gene['prod'])
			start = gene['pos'].split(' ')[0].split('..')[0]
			end = gene['pos'].split(' ')[0].split('..')[1]
			strand = gstrand(gene['pos'].split('[')[1].split(']')[0])
			product = gene['prod']

			if product == '':
				product = 'hypothetical protein'

			anno = {
					'start' : start,
					'end' : end,
					'strand' : strand,
					'locus_tag' : 'M'+ str(i),
					'type' : 'others',
					'color' : color,
					'description' : 'Location: '+gene['pos'].split(' ')[0]+' ('+gene['pos'].split(' ')[2]+' bp)<br>Type: ' +name +'<br>Detail: '+ product
				}
			if strand == 1:
				mapzlist.append(anno)
			else:
				mapflist.append(anno)
			i += 1

		head = 'var borders = [];\nvar tta_codons = [];\nvar orfs ='
		s = genelist[0]['pos'].split(' ')[0].split('..')[0]
		e = genelist[-1]['pos'].split(' ')[0].split('..')[1]

		gcdict = calculate_gc(fasta_file, int(s), int(e), 500, 50)
		with open(gcmap, 'r') as original_file:
			original_content = original_file.read()
		with open(gcjson,'w') as gein2:
			gein2t = 'var jsonData = ' + str(gcdict)+';'
			gein2.write(gein2t)
			gein2.write(original_content) 

#		with open(gcjson,'w') as gein2:
#			json.dump(gcdict, gein2, indent=4)

		maps = str(mapzlist)+';\nvar orfs2 ='+str(mapflist)+';\nvar clusterf2 = { start: '+s+', end: '+ \
				  e+', idx: 1, orfs: orfs, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nvar clusterr2 = { start: '+ s+', end: '+ \
				  e+', idx: 2, orfs: orfs2, borders: borders, tta_codons:tta_codons,\
				  label: \'\', unordered: true };\nsvgene.drawClusters("'+regijs+'", [clusterf2, clusterr2], 50, 920);'
		with open(mapfile,'w') as map_file:
			map_file.write(head + maps)

		with open(viewfile, 'r') as file:
			file_content = file.read()
		new_content = file_content.replace('XXXX', regijs)
		with open(htmlfile, 'w') as file:
			file.write(new_content)

	return ICEss

def copy_files(source_dir, destination_dir):
    if os.path.isfile(source_dir):
        shutil.copy(source_dir, destination_dir)
    elif os.path.isdir(source_dir):
        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)
        for item in os.listdir(source_dir):
            source_item = os.path.join(source_dir, item)
            destination_item = os.path.join(destination_dir, item)
            if os.path.isdir(source_item):
                copy_files(source_item, destination_item)
            else:
                shutil.copy2(source_item, destination_item)

def delete_folders_starting_with_keyword(dir, keyword):
	for dirpath, dirnames, filenames in os.walk(dir, topdown=False):
		for dirname in dirnames:
			if dirname.startswith(keyword):
				folder_to_remove = os.path.join(dirpath, dirname)
				shutil.rmtree(folder_to_remove)

def getfasta(runID,resultdir,id_dict,key,s,e,stag,etag):

	fafile = os.path.join(tmp_dir, runID, runID+'.fa')	
	faafile = os.path.join(tmp_dir, 'gbk', runID+'.faa')

	outfa = os.path.join(resultdir,key+'.fa')
	outfaa = os.path.join(resultdir,key+'.faa')

	seq_record = SeqIO.read(fafile, "fasta")
	with open(outfa, "w") as output_handle1:
		sequence = seq_record.seq[int(s)-1:int(e)]
		ID = '_'.join(seq_record.id.split('_')[-2:])
		seq_record.description = ''
		seq_record.seq = sequence
		seq_record.id = id_dict[ID] + ' ' + s +'-'+e
		SeqIO.write(seq_record, output_handle1, "fasta")

	faa_records = SeqIO.parse(faafile, "fasta")
	with open(outfaa, "w") as output_handle2:
		for faa_record in faa_records:
			seq_id = getnum(faa_record.id)
			if int(stag) <= seq_id <= int(etag):
				SeqIO.write(faa_record, output_handle2, "fasta")

def _meta(runID,infile):

	resultdir = os.path.join(workdir, 'result', runID)
	if not os.path.exists(resultdir):
		os.system('mkdir '+resultdir)
	jsback = os.path.join(workdir,'script','js')

	id_dict = rename(runID,infile)
	chosenfa = prescan(runID)
	drawout,spdict,report = Taxonomy(runID)
#	spdict = {}
	basefile = getbase(runID)
#	copy_files(drawout, resultdir)
	copy_files(basefile, resultdir)
	copy_files(report, resultdir)

	ICEsum = os.path.join(tmp_dir,runID, runID+'_ICEsum.json')
	newIDfa = os.path.join(tmp_dir, runID, runID+'_newID.fa')
	i = 1 
	ICEsumlist = []
	for seq_record in SeqIO.parse(newIDfa, "fasta"):
		if seq_record.id in chosenfa:
			sprunID = runID + '_' + seq_record.id
			seqfa = str(seq_record.seq)
			newfolder = os.path.join(tmp_dir, sprunID)			
			if not os.path.exists(newfolder):
				os.makedirs(newfolder)
			spfa = os.path.join(newfolder, sprunID+'.fa')
			with open(spfa,'w') as outfa:
				outfa.write(">%s\n%s\n" % (sprunID,seqfa))

			final_dir = os.path.join(tmp_dir,sprunID,'result')
			if not os.path.exists(final_dir):
				os.makedirs(final_dir)

			prokkanno(sprunID,spfa)
			ICEss = get_map(sprunID,spdict,id_dict)
			copy_files(final_dir, resultdir)

			if ICEss:
				for key,value in ICEss.items():
					[s,e,stag,etag] = value.split('|')
					lengt = int(e) - int(s) + 1
					ICEs = {
						'id' : str(i),
				        'seqid': id_dict[seq_record.id],
				        'species':spdict[seq_record.id],
#						'species':'',
				        'location': s+'..'+e,
				        'length': lengt,
				        'detail': key
				    }
					ICEsumlist.append(ICEs)
					getfasta(sprunID,resultdir,id_dict,key,s,e,stag,etag)
					i += 1 

	with open(ICEsum,'w') as ice_file:
		json.dump(ICEsumlist, ice_file, indent=4)

	copy_files(ICEsum, resultdir)
	jsdir = os.path.join(resultdir,'js')
	copy_files(jsback, jsdir)
	delete_folders_starting_with_keyword(tmp_dir, runID)
