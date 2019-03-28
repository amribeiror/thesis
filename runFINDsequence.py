import os
from os import listdir
import pandas as pd
from Bio import SeqIO
from linecache import getline
import sys

output = sys.argv[1]
input = sys.argv[2]

out_fastaL22 = output + "/fasta/LSU_22_5modelHMM.fa"
out_fastaL23 = output + "/fasta/LSU_23_jhmmer.fa"
out_fastaL24 = output + "/fasta/LSU_24_5modelHMM.fa"
out_xlsx= output + "/excel/ProOF.xlsx"
out_xlsx2= output + "/excel/plusANDminus.xlsx"
out_fafile= output + "/fasta/ribo_prot_nohitsonjhmmer.fa"

def HMM_to_SEQ():
    #in_directory = "../logs/hmmsearch/tblout"
    in_directory = input + "/proteomes"
    hit_directory = output + "/hmmsearch"
    out_xls = output + "/excel/hmmsearch.xlsx"
    out_fafile = output + "/fasta/hmmsearch_newnewnew_nohits.fa"
    out_nohits = output + "/excel/nohitscaet.xlsx"
    s = set()
    dL22 = {}
    dL23 = {}
    dL24 = {}
    unk = set()

    lista = []

    df = pd.DataFrame(columns=['QUERY', 'ORGANISM', 'TARGET', 'E-VALUE', 'TARGET_SEQUENCE', 'DESCRIPTION'])
    dg = pd.DataFrame(columns=['QUERY', 'ORGANISM'])

    for proteome in sorted(os.listdir(in_directory)) :
        org = proteome.split('_', 1)[0]
        with open(os.path.join(in_directory, proteome), 'r') as f :
            for record in SeqIO.parse(f, 'fasta'):
                for hit in os.listdir(hit_directory):
                    nohit = 'false'
                    if hit.startswith(org):
                        query = hit.split('_', 1)[1].strip('.txt')

                        if query + ':' + org + ':true' not in lista and query + ':' + org + ':false' not in lista :
                            lista.append(query + ':' + org + ':false')

                        with open(os.path.join(hit_directory, hit), 'r') as g :
                            #g.seek(-1, 2)  # go to the file end.
                            #eof = g.tell()  # get the end of file location
                            #g.seek(0, 0)  # go back to file beginning

                            #while g.tell() != eof:

                            #assert n > 0

                            #while True:

                                #chunk = g.read(n)
                            for i, l in enumerate(g):


                                #print(l, 'this is a line')
                                count = i + 1
                                if count == 4:
                                    if l.startswith('#'):
                                        #print('No homologue found in ' + org + ' for ' + query)
                                        if org + '_' + query not in unk :
                                            unk.add(org + '_' + query)
                                            dg = dg.append({"QUERY": query, "ORGANISM": org}, ignore_index=True)

                                #print(proteome, 'this is the proteome')

                                #for line in g :
                                if not l.startswith('#') :

                                    cols = [x for x in l.strip().split(' ') if x]
                                    if cols[0] == record.name :
                                        #print(record.name, 'this is record name')
                                        #print(cols[0], 'this is target')
                                        #nohit = '0'
                                        if str(record.seq) not in s :
                                            if len(cols) > 19:
                                                cols[18] = ' '.join(cols[18:])
                                                #    if it's < 19, we have no description columns, so use an empty string
                                                # instead
                                            elif len(cols) < 19:
                                                cols.append('')
                                                assert len(cols) == 19

                                                # # assign parsed column data into qresult, hit, and hsp dicts
                                                # qresult = {}
                                                # qresult['id'] = cols[2]  # query name
                                                # qresult['accession'] = cols[3]  # query accession
                                                # hit = {}
                                                # hit['id'] = cols[0]  # target name
                                                # hit['accession'] = cols[1]  # target accession
                                                # hit['evalue'] = float(cols[4])  # evalue (full sequence)
                                                # hit['bitscore'] = float(cols[5])  # score (full sequence)
                                                # hit['bias'] = float(cols[6])  # bias (full sequence)
                                                # hit['domain_exp_num'] = float(cols[10])  # exp
                                                # hit['region_num'] = int(cols[11])  # reg
                                                # hit['cluster_num'] = int(cols[12])  # clu
                                                # hit['overlap_num'] = int(cols[13])  # ov
                                                # hit['env_num'] = int(cols[14])  # env
                                                # hit['domain_obs_num'] = int(cols[15])  # dom
                                                # hit['domain_reported_num'] = int(cols[16])  # rep
                                                # hit['domain_included_num'] = int(cols[17])  # inc
                                                # hit['description'] = cols[18]  # description of target
                                                # #print(hit['description'])


                                                #print(record.id, 'this is record id')

                                            s.add(str(record.seq))
                                            lista = [w.replace(query + ':' + org + ':false', query + ':' + org + ':true') for w in lista]
                                            #lista.replace('false', 'true')
                                            print('A putative homologue of ' + query + ' was found in ' + org)
                                                #with open(out_file, 'a') as h :
                                                #    h.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (query, org, record.name, cols[4], str(record.seq), cols[18]))
                                            df=df.append({"QUERY": query, "ORGANISM": org, "TARGET": record.name, "E-VALUE": cols[4], "TARGET_SEQUENCE": str(record.seq), "DESCRIPTION" : cols[18]}, ignore_index=True)
                                            if query == 'uL22m':
                                                dL22[record.name] = str(record.seq)
                                            elif query == 'uL23m':
                                                dL23[record.name] = str(record.seq)
                                            if query == 'uL24m':
                                                dL24[record.name] = str(record.seq)
                                #print("file finished")
                    #if "" == l:
                    #    print("file finished")
    #print(lista)
    with open(out_fafile, 'w') as j:
        for nohits in lista :
            if 'false' in nohits:
                j.write('No putative homologues were found for:%s\n' % (nohits))

                                #         #print("file finished")

                                #         with open(out_fafile, 'a') as j:
                                #             j.write('>%s\n%s\n' % (query, proteome))
                                # #break


    df.to_excel(out_xls)
    df.to_excel(out_nohits)
    with open(out_fastaL22, 'w') as L22:
        for k, v in dL22.items():
            L22.write('>%s\n%s\n' % (k, v))
    with open(out_fastaL23, 'w') as L23:
        for k, v in dL23.items():
            L23.write('>%s\n%s\n' % (k, v))
    with open(out_fastaL24, 'w') as L24:
        for k, v in dL24.items():
            L24.write('>%s\n%s\n' % (k, v))

HMM_to_SEQ()