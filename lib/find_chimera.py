#! /bin/env python

import sys
sys.path.append("./")
from get_corename import get_corename

##########################################################################
def set_find_chimera_args():
    find_chimera_args={}
    find_chimera_args['num_of_units']=50
    find_chimera_args['coverage'] = {'min':0.2,'max':0.7}
    find_chimera_args['aligned_length'] = {'min':150}
    return find_chimera_args
  

def generate_chimera_info_all_posi(num_of_units,gene_length):
    all_posi_array=[]
    length_per_unit_floated = gene_length/float(num_of_units)
    for i in range(1,num_of_units+1):
        all_posi_array.append(int(length_per_unit_floated*i))
    return(all_posi_array,length_per_unit_floated)

 
def get_chimera_info(query,subject,start,end,evalue,chimera_info_gene_dict,gene_length,find_chimera_args,
                    is_no_corename):
    # query, subject do not represent the real query and subject and could be switched
    if not 'evalue' in chimera_info_gene_dict:
        chimera_info_gene_dict['evalue']={}
        chimera_info_gene_dict['hit']={}
    for i in chimera_info_gene_dict['all_posi']:
        if start <= i < end:
            chimera_info_gene_dict=fill_in_chimera_info_gene_dict(chimera_info_gene_dict,evalue,subject,i,is_no_corename)
        elif i >= end:
            break
    return(chimera_info_gene_dict)


def fill_in_chimera_info_gene_dict(chimera_info_gene_dict,evalue,subject,i,is_no_corename):
    subject_corename = get_corename(subject,is_no_corename)
    if not i in chimera_info_gene_dict['evalue']:
        chimera_info_gene_dict['evalue'][i] = {}
    chimera_info_gene_dict['evalue'][i][subject_corename] = evalue
    #chimera_info_gene_dict['hit'][i] = subject_corename
    '''    
    if not subject_corename in chimera_info_gene_dict['hit']:
        chimera_info_gene_dict['hit'][subject_corename]=[]
    chimera_info_gene_dict['hit'][subject_corename].append(i)
    '''
    return(chimera_info_gene_dict)


def parse_chimera_info(chimera_info,outfile):
    if not chimera_info:
        return # if no chimera_info is generated, return
    fh=open(outfile,'w')
    for gene in chimera_info.keys():
        if (not 'evalue' in chimera_info[gene]) or (not chimera_info[gene]['evalue']):
            continue
        
        chimera_info[gene]['hit'] = {}
        for posi in sorted(chimera_info[gene]['evalue']):
            a = sorted(chimera_info[gene]['evalue'][posi].items(), key=lambda d: d[1])
            chimera_info[gene]['hit'][posi] = a[0][0]
                
        effective_hits={}
        for posi,hit in chimera_info[gene]['hit'].iteritems():
            if not hit in effective_hits:
                effective_hits[hit]=[]
            effective_hits[hit].append(posi)

        num_of_effective_hits = 0
        for posi in effective_hits[hit]:
            values = chimera_info[gene]['hit'].values()
            best_hits = filter (lambda x: values.count(x)>=10, values)

        for hit in dict(zip(best_hits,['']*len(best_hits))):
            counter_4_an_effective_hit=0
            other_hits = filter (lambda x: x != hit, best_hits)
            #other_hits = filter (lambda x: x != hit, effective_hits.keys())
            #for posi in effective_hits[hit]:
            for posi in effective_hits[hit]:
                OtherHitsIn_array = []
                OtherHitsIn_array = filter(lambda hit: hit in chimera_info[gene]['evalue'][posi], other_hits)
                if not OtherHitsIn_array:
                    counter_4_an_effective_hit += 1
            if (counter_4_an_effective_hit >= 10) or (counter_4_an_effective_hit * chimera_info[gene]['length_per_unit_floated'] >= 150):
                num_of_effective_hits += 1
    
        if num_of_effective_hits >= 2:
            fh.write(gene+'\n')
            fh.write('\t'.join(effective_hits)+'\n')
            print chimera_info[gene]['hit']
            #print chimera_info[gene]['evalue']
    fh.close()


