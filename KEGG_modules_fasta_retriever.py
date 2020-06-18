#!/usr/bin/env python
# coding: utf-8

# In[163]:


#!/usr/bin/env python


# In[1]:


import requests


# In[2]:


import json
import re
import functools
import operator
import os


# In[3]:


from collections import defaultdict


# In[5]:


#Creating organism list with Euk/Pro info from KEGG


# In[6]:


organisms=requests.get("http://rest.kegg.jp/list/organism/")


# In[7]:


organisms_data=(organisms.content).decode("utf-8")


# In[8]:


org_list=organisms_data.split("\n")


# In[9]:


org_dict={}
for line in org_list[0:-1]:
    org_dict[line.split("\t")[1]] = [str(line.split("\t")[3].split(";")[0])]


# In[10]:


org_dict_sub={}
for line in org_list[0:-1]:
    org_dict_sub[line.split("\t")[1]] = [str(line.split("\t")[3].split(";")[1])]


# In[11]:


#Retrieving all the module IDs from KEGG


# In[ ]:


#ko00002.keg file downloaded from https://www.genome.jp/kegg-bin/get_htext
#contains all modules of KEGG database


# In[12]:


mod_data=open("ko00002.keg", "r")


# In[13]:


mod_lines=mod_data.readlines()


# In[14]:


mod_exp="M\d{5}"
module_names_list=[]
for line in mod_lines:
    if len(line.split())>2:
        mod_id=line.split()[1]
        if re.search(mod_exp, mod_id):
            module_names_list.append(mod_id)
        


# In[15]:


#Function module to ortholog
#input: module id
#output: orthologs list from that module


# In[16]:


def module_to_ortholog(module_id):
    mod_file_1=requests.get("http://rest.kegg.jp/get/%s" %str(module_id))
    mod_file=(mod_file_1.content).decode("utf-8").split()
    ort_exp="K\d{5}"
    mod_file2=[]
    x=False
    for i in range(len(mod_file)-1):
        if mod_file[i]=="ORTHOLOGY":
            x=i
    for j in range(x, len(mod_file)-1):
        if re.findall(ort_exp, mod_file[j]):
            y=re.findall(ort_exp, mod_file[j])
            mod_file2.append(y)
    ortholog_list=functools.reduce(operator.iconcat, mod_file2, [])
    return ortholog_list


# In[17]:


#Function ortholog to gene
#input: orthology id 
#output: genes list (with specie abbrevation)


# In[18]:


def ortholog_to_gene(ortholog_id):
    ort_file=requests.get("http://togows.org/entry/kegg-orthology/%s/genes.json" % str(ortholog_id))
    ort_genes1=ort_file.json()
    ort_genes=ort_genes1[0]
    ort_genes_list=[]
    for key in ort_genes:
        ort_genes_list.append([key])
    return ort_genes_list


# In[ ]:


#Function ortholog to gene
#input: orthology id 
#output: genes list (with specie abbrevation and gene_id)


# In[19]:


def ortholog_to_gene_id(ortholog_id):
    ort_file=requests.get("http://togows.org/entry/kegg-orthology/%s/genes.json" % str(ortholog_id))
    ort_genes1=ort_file.json()
    ort_genes=ort_genes1[0]
    ort_genes_list=[]
    for key in ort_genes:
        for item in ort_genes[key]:
            ort_genes_list.append(([key], item))
    return ort_genes_list


# In[20]:


#Function to tell Eukaryote/Prokaryote from the species abbreviation of gene code
#input: gene specie abbrevation
#output: ["Eukaryote"] / ["Prokaryote"]


# In[21]:


def tax_from_gene(gene_id):
        gene_id_2=str(gene_id)[2:-2]
        tax=org_dict[gene_id_2]
        return tax


# In[ ]:


#Function to retrieve sub taxonomy of the gene (if necessary)


# In[22]:


def sub_tax_from_gene(gene_id):
    tax=org_dict_sub[gene_id]
    return tax


# In[ ]:


#Creating a subset of list of modules that has only prokaryotic genes


# In[23]:


final_modules=[]
tax_list=[]
for line in module_names_list:
    orthologs=module_to_ortholog(line)
    mod_condition=False
    for ortholog in orthologs:
        genes=ortholog_to_gene(ortholog)
        if mod_condition==True:
            break
        try:
            for gene in genes:
                taxon=tax_from_gene(gene)
                if taxon==["Eukaryotes"]:
                    mod_condition=True
                    break
                elif taxon==["Prokaryotes"]:
                    continue
                else:
                    print(tax_from_gene(gene)+"error")
        except KeyError:
            continue
    if mod_condition==False:
        final_modules.append(line)


# In[ ]:


#Retrieving fasta files from KEGG database


# In[77]:


for line in final_modules:
    orthologs=module_to_ortholog(line)
    for ortholog in orthologs:
        genes=ortholog_to_gene_id(ortholog)
        ortholog_file=open("%s_%s.fas" %(str(line)[0:-1], str(ortholog)), "w")
        for gene in genes:

            file=requests.get("http://rest.kegg.jp/get/%s:%s/aaseq" %(str(gene[0])[2:-2],str(gene[1])))
            fasta_1=(file.content).decode("utf-8")
            ortholog_file.write(fasta_1)
    ortholog_file.close()

