#!/usr/bin/python
import argparse
import os
import re
import pdb
#from guppy import hpy
#import resource
#from profilestats import profile
#import cProfile
#import line_profiler
#import time


def vcf2json_alt(v, vcf_content ):
	vf = open(v, 'r')
	status_bar = ''

	for aline in vf:
		if(not re.match('^#', aline)):
			aline = re.sub('\s+$', '', aline)
			content = aline.split('\t')
			aVar = {}
			#full_loc = conv[content[0]]
			#alt_loc = conv[content[0]].split('_')
			var_chr = content[0]
			
			var_p   = content[1]
			#var_rsid  = content[2]
			var_ref = content[3]
			var_alt = content[4]
			if( var_chr == 'chrX' or var_chr == 'X'):
				var_chr = 'chr23'
			elif (var_chr =='chrY' or var_chr == 'Y'):
				var_chr = 'chr24'
			elif(var_chr =='MT'or var_chr =='M'):
				var_chr = 'chr25'
			if status_bar != var_chr:
				status_bar = var_chr
				print(status_bar)
				#print(status_bar, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, ' MB')
			var_id = var_chr+':'+var_p
			#alt_allele = re.split(',', var_alt)
			#num_allele = len(re.split(',', var_alt))
			#grp_size = []
			total_size = 0
			#cln_allele = {}
			info = content[-1].split(';')
			
			tmp = {}
			for ele in info:  ## parse the info column of VCF file
				pair = ele.split('=')
				
				if(len(pair) ==2):
					if pair[0] == 'CAF':
						#all_freq = re.findall('\d{1}|\d\.\d+', pair[1])
						all_freq = pair[1].split(',')
						all_freq = [float(x) for x in all_freq if x != '.']
						refMaf = all_freq[0]
						all_freq = sorted(all_freq)
						if len(all_freq) >=2 and refMaf == all_freq[-2]:
							tmp['RMA'] = '1'
						if len(all_freq) >=2:
							maf = all_freq[-2]
							maf = round(maf, 3)
							str_maf = str(maf)
							if maf < 0.001:
								str_maf = '0.0'
							#else:
							#	str_maf = str(round(maf, 3))
							
							tmp[pair[0]] = str_maf
					else:
						tmp[pair[0]] = pair[1]
			tmp['o'] = var_alt
			tmp['r'] = var_ref
			
			if var_id not in vcf_content:
				vcf_content[var_id] = [tmp]
			else:
				vcf_content[var_id].append(tmp)
	vf.close()
	return('chr'+status_bar)

def inUCSC(fl, ucsc):
	uf = open(fl, 'r')
	header = uf.readline()
	for u in uf:
		aline = u.strip().split('\t')
		ucsc[aline[3]] = 1

def testdbSNP(infle, ucsc):

	oridata = {}
	vcf_small_file = vcf2json_alt(infle,  oridata)
	out_f = open(vcf_small_file.strip()+'_test.json', 'w')

	#meta = []
	for var,all_info in oridata.items():
		id_inf = ''
		id_json = ''
		#var_pos = ''
		delim = re.compile('[:_]')
		id_e = delim.split(var)
		id_e[0] = re.sub('chr', '',id_e[0])
		var_pos = id_e[-1]
		id_json = '{\"_id\":{\"c\":'+id_e[0]+',\"p\":'+id_e[1]+'},'
		json_line = []
		for info in all_info:
			meta = []
			meta.append('\"i\":'+'\"rs'+info['RS'] +'\"')
			if 'CAF' in info:
				meta.append('\"maf\":'+info['CAF'])
			alt_allele = info['o'].split(',')
			f_o = []
			for a in alt_allele:
				f_o.append('\"'+a+'\"')
			meta.append('\"o\":[' + ','.join(f_o) +']')
			if (info['RSPOS'] != var_pos and info['VC'] != 'DIV') or (info['VC'] == 'DIV' and int(info['RSPOS']) != (int(var_pos)+1)):
				meta.append('\"opos\":'+ info['RSPOS'])
			meta.append('\"r\":'+'\"'+info['r'] +'\"')
			if 'RMA' in info:
				meta.append('\"rma\":'+info['RMA'])
			meta.append('\"s\":'+info['SAO'])
			info['VC'] = re.sub('DIV', 'INDEL', info['VC'])
			meta.append('\"t\":['+'\"'+info['VC'] +'\"]')
			if 'rs'+info['RS'] in ucsc:
			#if info['ucsc'] ==1:
				meta.append('\"ucsc\":1')

			#json_line = json_line + ',{' + ','.join(meta)+ '}'
			json_line.append('{' + ','.join(meta) + '}')
		out_line = id_json +'\"f\":[' + ','.join(json_line) + ']}'
		out_f.write(out_line+'\n')

	out_f.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input')
	args = parser.parse_args()
	infls = open(args.input, 'r')
	
	ucsc = {}
	inUCSC('GRCh38_ucsc_commonSNPs144.bed', ucsc)

	for f in infls:
		testdbSNP(f.strip(), ucsc)
	infls.close()
	ucsc = {}

