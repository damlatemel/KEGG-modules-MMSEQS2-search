#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python


# In[19]:


import csv


# In[20]:


import re


# In[21]:


import sys


# In[23]:


def fasta_reader(filename):    
    fas_file=open(filename, "r")
    seqs={}
    seq_fragments=[]
    for line in fas_file:
        if line.startswith(">"):
            seq_id=line.rstrip()[1:]
            if seq_fragments:
                sequence=''.join(seq_fragments)
                seqs[seq_id]=sequence
            seq_fragments=[]
        else:
            seq=line.rstrip()
            seq_fragments.append(seq)
    if seq_fragments:
        sequence= ''.join(seq_fragments)
        seqs[seq_id]=sequence
    fas_file.close()
    return seqs


# In[24]:


def noclass_contig(file_name_arg):
    set_tsv=open(file_name_arg)
    set_file=csv.reader(set_tsv, delimiter="\t")
    euk_exp="Eukaryot"
    bac_exp="Bacteria"
    arc_exp="Archae"
    vir_exp="Virus"
    contig_dict={}
    for row in set_file:
        if re.search(euk_exp, str(row)):
            contig_dict[str(row[0])]="Eukaryote"
        elif re.search(bac_exp, str(row)):
            continue
            
        elif re.search(arc_exp, str(row)):
            continue
            
        elif re.search(vir_exp, str(row)):
            continue
            
        else:
            contig_dict[str(row[0])]="Unknown"
            
    return contig_dict


# In[38]:


def fasta_subseter(fasta_file_name, contig_dict, taxon, output_file):
    fasta_file=open(output_file, "w")
    fasta_list=[]
    for line in fasta_file_name:
        tax_id=str(line).split("|")[1]
        if str(tax_id) in contig_dict and contig_dict[str(tax_id)]==str(taxon):
            fasta_file.writelines("\n>%s \n%s" % (str(line), str(fasta_file_name[line])))
    fasta_file.close()
    return fasta_file
        


# In[ ]:


fasta_file_name=sys.argv[1]


# In[ ]:


contig_name=sys.argv[2]


# In[ ]:


taxon_group=sys.argv[3]


# In[ ]:


output_file_name=sys.argv[4]


# In[ ]:


fasta_files=fasta_reader(fasta_file_name)


# In[ ]:


contigs_dict=noclass_contig(contig_name)


# In[39]:


fasta_subseter(fasta_files,contigs_dict, taxon_group, output_file_name )


# In[ ]:




