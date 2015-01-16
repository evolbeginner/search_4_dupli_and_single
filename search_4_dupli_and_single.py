#! /usr/local/bin/python

# To identify duplicates, singletons and/or duplciated gene pairs based on reciprocal BLAST hits.
# Author: Sishuo Wang from the department of Botany, the University of British Columbia
# Email: sishuowang@hotmail.ca
# Note: some packages need to be installed in advance.

#######################################################################
import sys
import getopt
import re
import os
from Bio import SeqIO
import operator
import time
# import selfbasics

# could be ignored if the option '--find_chimera' is not specified
sys.path.append(os.path.join(os.path.dirname(__file__),'lib'))
import find_chimera
from get_corename import get_corename


#######################################################################
duplicates={}
singletons={}
non_singletons={}
pair_with={}
reciprocal_pairs={}
pair_list_excluded=[]
pairs_excluded={}
seq_objs={}
resort_swi=''
upstreamPointOnS4Q={}
find_chimera_args={}
chimera_info={}

#######################################################################
def read_param():
    outdir=''
    is_force=False
    duplicate_evalue=float(1e-10)
    singleton_evalue=float(1e-3)
    min_identity=0
    min_bit_score=0
    min_coverage=0
    coverage_query_swi=0
    seq_file_format='fasta'
    is_duplicate_pairs=False
    best_evalue_items=[]
    no_TD_swi=''
    resort_swi=''
    is_complement=''
    is_no_corename=False
    is_find_chimera=False
    find_chimera_args={}
    alter_min_AlignedLength=200
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help","seq_file=","seq_file_format=","blast_file=","outdir=","force","duplicate_evalue=","singleton_evalue=","min_identity=","min_bit_score=","min_coverage=","alter_min_AlignedLength=","coverage_query", "duplicate_pairs", "pair_list_excluded=", "best_evalue_item=", "no_TD", "resort", "complement", "no_corename", "find_chimera"])
    except getopt.GetoptError:
        print "Illegal params!"
        show_help()
    for op, value in opts:
        if op == '--seq_file':
            seq_file=os.path.expanduser(value)
        elif op == '--seq_file_format':
            seq_file_format=value
        elif op == '--blast_file':
            blast_file=os.path.expanduser(value)
        elif op == "--outdir":
            outdir=os.path.expanduser(value)
        elif op == "--force":
            is_force=True
        elif op == "--duplicate_evalue":
                duplicate_evalue=float(value)
        elif op == "--singleton_evalue":
            singleton_evalue=float(value)
        elif op == "--min_bit_score":
            min_bit_score=float(value)
        elif op == "--min_identity":
            min_identity=float(value)
        elif op == "--min_coverage":
            min_coverage=float(value)
        elif op == "--alter_min_AlignedLength":
            alter_min_AlignedLength=int(value)
        elif op == "--coverage_query":
            coverage_query_swi=1
        elif op == "--duplicate_pairs":
            is_duplicate_pairs=True
        elif op == "--pair_list_excluded":
            pair_list_excluded.append(os.path.expanduser(value))
        elif op == "--best_evalue_item":
            for i in value.split(','):
                best_evalue_items.append(i)
        elif op == "--no_TD":
            no_TD_swi=True
        elif op == "--resort":
            resort_swi = 1
        elif op == "--complement":
            is_complement = 1
        elif op == "--no_corename":
            is_no_corename = True
        elif op == "--find_chimera":
            is_find_chimera = True
        elif re.search('^--?h(elp)?$',op):
            show_help()
#    ================================    #
    if not sys.argv[1:]:
        show_help()
    if not outdir:
        print "outdir has to be specified! Exiting ......"
    else:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if os.path.exists(outdir) and is_force:
            rm_dir_cmd = "rm -rf %s" % outdir
            os.system(rm_dir_cmd)
            os.makedirs(outdir)

    if is_find_chimera:
        find_chimera_args = find_chimera.set_find_chimera_args()

    out_param=open(outdir+'/Params', 'w')
    timeArray = time.localtime(time.time())
    otherStyleTime = time.strftime("%Y-%m-%d %H:%M:%S", timeArray)
    out_param.write(otherStyleTime+'\n')
    for i in sys.argv[1:]:
        out_param.write(i+'\n')
    out_param.close()
#    ================================    #
    return (blast_file, seq_file, seq_file_format, outdir, duplicate_evalue, singleton_evalue, min_identity, min_bit_score, min_coverage, alter_min_AlignedLength, is_duplicate_pairs, pair_list_excluded, coverage_query_swi, best_evalue_items, no_TD_swi, resort_swi, is_complement, is_no_corename, find_chimera_args)


def read_seq_file(seq_file, seq_file_format='fasta'):
    class Seq:
        def __init__(self):
            self.length = len(seq_record.seq)
    print "reading sequence file"
    handle = open(seq_file, "r")
    for seq_record in SeqIO.parse(handle, seq_file_format):
        seq_objs[seq_record.id] = Seq()
    handle.close()


def read_blast_file(blast_file, duplicate_evalue, singleton_evalue, min_identity, min_bit_score, min_coverage, alter_min_AlignedLength, coverage_query_swi, best_evalue_items, find_chimera_args={}):
    highest_bit_score={}
    blast_file_handle=open(blast_file,'r')
    chimera_info={}
    print "reading BLAST file"

    for line in blast_file_handle:
        line = line.rstrip('\n\r')
        line2 = line.split('\t')
        [query,subject,identity,aligned_length,query_start,query_end,subject_start,subject_end,evalue,bit_score] = operator.itemgetter(0,1,2,3,6,7,8,9,10,11)(line2)
        identity=float(identity)
        evalue=float(evalue)
        query_start,query_end=int(query_start),int(query_end)
        subject_start,subject_end=int(subject_start),int(subject_end)
        aligned_length=float(aligned_length) # must be float since in will be involved in ChuFa
        bit_score=float(bit_score)
        query_corename = get_corename(query,is_no_corename); subject_corename = get_corename(subject,is_no_corename)
        if (query_corename != subject_corename):
            if (evalue <= singleton_evalue):
                non_singletons[query_corename]=1
            if (evalue > duplicate_evalue):
                continue
            if (min_bit_score>0):  # whether 'min_bit_score' has been defined
                if (bit_score < min_bit_score):
                    continue
            if min_identity > 0:
                if identity < min_identity:
                    continue
            if min_coverage > 0: # defalut: 0
                if coverage_query_swi:
                    if (aligned_length/seq_objs[query].length) < min_coverage \
                        and aligned_length <= alter_min_AlignedLength:
                        continue
                else:
                    if aligned_length/min(seq_objs[query].length, seq_objs[subject].length) < min_coverage:
                        continue

            if find_chimera_args:
                #for index,gene in enumerate([query,subject]):
                for index,gene in enumerate([query]):
                    if not gene in seq_objs:
                        break
                    else:
                        if not gene in chimera_info:
                            chimera_info[gene]={}
                        if find_chimera_args['coverage']['min'] <= aligned_length/seq_objs[gene].length or \
                           find_chimera_args['aligned_length']['min'] <= aligned_length :
                            if not 'all_posi' in chimera_info[gene]:
                                chimera_info[gene]['all_posi'],chimera_info[gene]['length_per_unit_floated'] = find_chimera.generate_chimera_info_all_posi(find_chimera_args['num_of_units'],seq_objs[gene].length)
                            if index == 0:
                                chimera_info[gene] = find_chimera.get_chimera_info(gene,subject,query_start,query_end,evalue,chimera_info[gene],seq_objs[gene].length,find_chimera_args,is_no_corename)
                            elif index == 1:
                                chimera_info[gene] = find_chimera.get_chimera_info(gene,query,subject_start,subject_end,evalue,chimera_info[gene],seq_objs[gene].length,find_chimera_args,is_no_corename)
                            else:
                                raise MyError('index out of range')
                        
            if best_evalue_items:
                is_continue=0
                q_s = {'subject':subject,'query':query}
                for i in best_evalue_items:
                    if not q_s[i] in highest_bit_score:
                        highest_bit_score[q_s[i]]=bit_score
                    else:
                        if bit_score <= highest_bit_score[q_s[i]]:
                            is_continue=1
                if is_continue == 1:
                    continue

            duplicates[query_corename] = 1
            pair = [query_corename, subject_corename]
            for i in [0,1]:
                gene1 = pair[i]
                gene2 = pair[abs(1-i)]
                if not gene1 in pair_with:
                    pair_with[gene1]={}
                pair_with[gene1][gene2] = evalue
            if best_evalue_items:
                upstreamPointOnS4Q[query_corename]={}
                upstreamPointOnS4Q[query_corename]['start']=query_start
                upstreamPointOnS4Q[query_corename]['hit_start']=subject_start
                upstreamPointOnS4Q[query_corename]['paralog']=subject_corename
                if subject_start<subject_end:
                    upstreamPointOnS4Q[query_corename]['direction']='+'
                else:
                    upstreamPointOnS4Q[query_corename]['direction']='-'
    return(chimera_info)


def select_singletons():
    for i in seq_objs.keys():
        m=re.search(r'^(\S+)',i)
        if m:
            corename=get_corename(m.group(1),is_no_corename)
            if is_complement:
                if (not corename in duplicates):
                    singletons[corename]=1
            else:
                if (not corename in non_singletons):
                    singletons[corename]=1


def get_pairs_excluded(pair_list_arr, seperator, resort_swi):
    pairs_excluded = {}
    for list_file in pair_list_arr:
        ff = open(list_file, 'r')
        for line in ff:
            line = line.rstrip('\r\n')
            line_arr = line.split(seperator)
            line_arr = line_arr[0:2] # =========================
            if resort_swi:
                pair = "|".join(sorted(line_arr))
            else:
                pair = "|".join(line_arr)
            pairs_excluded[pair] = 1
    return (pairs_excluded)


def output(duplicates, singletons, pairs, outdir, is_duplicate_pairs, upstreamPointOnS4Q):
    print "outputting results"
    out1=open(outdir+'/duplicates.list', 'w')
    for i in sorted(duplicates.keys()):
        out1.write(i+'\n')
    out1.close()

    out2=open(outdir+'/singletons.list', 'w')
    for i in sorted(singletons.keys()):
        out2.write(i+'\n')
    out2.close()

    if is_duplicate_pairs:
        out3=open(outdir+'/pairs.list', 'w')
        for i in sorted(pairs.keys()):
            i = re.sub('\|', '\t', i)
            out3.write(i+'\n')
        out3.close()

    if upstreamPointOnS4Q:
        out4=open(outdir+'/upstreamPointOnS4Q','w')
        for k,v in sorted(upstreamPointOnS4Q.iteritems(),key=lambda d:d[0]):
            out4.write(k+"\t"+str(v['paralog'])+"\t"+str(v['direction'])+"\t"+str(v['start'])+"\t"+str(v['hit_start'])+"\n")
        out4.close()


###################################################################################################
def is_TD(gene1, gene2):
    ref = re.compile('^AT([0-9MC])G(\d+)', re.I)
    m1=ref.search(gene1)
    m2=ref.search(gene2)
    if m1.group(1) == m2.group(1) and float(m1.group(2))-float(m2.group(2)) <= 1000:
        return (True)
    else:
        return (False)

def get_pair(pair_with, pairs_excluded, no_TD_swi):
    has_processed={}
    paired={} # paired[gene1] = gene2
    for gene in pair_with:
        lowest_evalue = 1000
        paired_with = ''
        for paralog in pair_with[gene]:
            if lowest_evalue >= pair_with[gene][paralog]:
                lowest_evalue = pair_with[gene][paralog]
                paired_with = paralog
                paired[gene] = paired_with
    # get reciprocal pairs
    for gene in paired:
        if gene in has_processed:
            continue
        tmp_paired = paired[gene]
        has_processed[tmp_paired] = 1
        if gene == paired[tmp_paired]:
            a = "|".join(sorted([gene, tmp_paired]))
            if a in pairs_excluded:
                continue
            if no_TD_swi:
                if is_TD(gene, paired[gene]):
                    continue
            reciprocal_pairs[a] = 1


def show_help():
    print '''    Usage:  python2.7 search_4_dupli_and_single.py <--seq_file=sequence_file> <--blast_file=blast_result> <--outdir=output_directory> [options]
    Options:
    --force               If outdir exists, remove it and create a new one with the same name.
    --seq_file_format=    the format of the sequence file
                          default: fasta
    --singleton_evalue=   the minimum evalue for singletons
                          default: 1e-3
    --duplicate_evalue=   the maximum evalue for duplicates
                          default: 1e-10
    --min_identity=       the minimum identity of the alignment
    --min_bit_score=      the minimum bit score
                          default: 0
                          e.g. 50
    --min_coverage=       the minimum value of aligned sequence length divided by the minimum sequence length between query and subject
                          default: 0
                          e.g. 0.5
    --alter_min_AlignedLength
                          If the coverage does not pass the threshold, alternatively, this value will be used to see if the aligne_length is smaller than the value given. If so, sequence will not be considered as duplicates.
    --coverage_query      If this argument is specified, the minimum coverage will be calculated as the aligned sequence length divided by the sequence length of the query.
                          Default: OFF
    --best_evalue_item=   Only the hit(s) with the lowest evalue will be chosen.
                          available options: query, subject
                          Default: OFF
    --duplicate_pairs     to identify duplicate pairs?
    --pair_list_excluded= pair_list file(s) where gene pairs would be excluded in final reciprocal gene pairs 
    --no_TD               no tandem duplicates based on the gene locus
    --resort              resort gene pairs in the given lists (--pair_list_excluded)
    --is_complement       select duplicates first and all others will be defined as singletons
    --no_corename         no corename conversion
    --find_chimera        find chimeric genes
                          Default: OFF
    -h|--h|--help         show usage

    Note:
    1. Blast output file should be in the format of format 8 (http://edwards.sdsu.edu/labsite/index.php/ramys/238-blast-output-8).
    2. Biopython module needs to be installed.
'''
    sys.exit(0)

#######################################################################
(blast_file, seq_file, seq_file_format, outdir, duplicate_evalue, singleton_evalue,
 min_identity, min_bit_score, min_coverage, alter_min_AlignedLength, is_duplicate_pairs, pair_list_excluded, coverage_query_swi,
 best_evalue_items, no_TD_swi, resort_swi, is_complement, is_no_corename, find_chimera_args) = read_param()

read_seq_file(seq_file, seq_file_format)

chimera_info = read_blast_file(blast_file, duplicate_evalue, singleton_evalue, min_identity, min_bit_score, min_coverage, alter_min_AlignedLength, coverage_query_swi, best_evalue_items, find_chimera_args)

select_singletons()

find_chimera.parse_chimera_info(chimera_info,os.path.join(outdir,'chimera_output'))

if is_duplicate_pairs:
    if pair_list_excluded:
        pairs_excluded = get_pairs_excluded(pair_list_excluded, "\t", resort_swi)
    get_pair(pair_with, pairs_excluded, no_TD_swi)

output(duplicates, singletons, reciprocal_pairs, outdir, is_duplicate_pairs, upstreamPointOnS4Q)


