#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#MMSEQS2 output filtering


# In[1]:


import pandas as pd


# In[2]:


pd.set_option('display.max_colwidth', None)


# In[3]:


#Make a dataframe out of the output file of MMSEQS2 search


# In[3]:


def output_to_dataframe(filename):
    eukbook=open(filename).readlines()
    hits=[]
    for line in eukbook:
        hit=dict()
        hit["contig_id"]=line.split("|")[1]
        hit["sample_id"]=line.split("|")[1].split("_")[0]
        hit["protein_id"]=line.split("|")[0]
        hit["module_id"]=line.split("\t")[-1][0:6]
        hit["ortholog_id"]=line.split("\t")[-1][7:13]
        hit["gene_id"]=line.split("\t")[1].split(":")[1]
        hit["E_value"]=float(line.split("\t")[-4])
        hit["original_line"]=line
        hits.append(hit)
    eukbook_df=pd.DataFrame(hits, columns=["sample_id", "contig_id","protein_id", "module_id", "ortholog_id", "gene_id","E_value", "original_line"])
    return eukbook_df


# In[9]:


#FUnction: Removing duplicates of hits that are in the
#1) same contig, same module and same ortholog (each ortholog has many genes inside which causes the repetition in the hits)
#2) same sample, same module and same ortholog


# In[10]:


def remove_duplicates(filename, sample_level):
    eukbook_filtered=filename.sort_values(by=["E_value"]).drop_duplicates(subset=[sample_level, "module_id", "ortholog_id"])
    return eukbook_filtered


# In[73]:


def df_grouper(filename, sample_level):
    #Retrieving the ortholog number in each module from previously created text file
    mod_count=open("mod_count.txt", "r").readlines()
    mod_dict={}
    for line in mod_count:
        mod_dict[line.split(",")[0][2:-1]]=int(line.split(",")[1])
    #Grouping the data frames
    df_grouped=filename.groupby([sample_level, "module_id"], as_index=False).agg({"ortholog_id":["count"], "protein_id":["nunique"]})
    df_grouped.columns=df_grouped.columns.droplevel(0)
    df_grouped.columns=[sample_level,'module_id', 'ortholog_number', 'protein_number']
    df_grouped=df_grouped.sort_values(by=["ortholog_number"], ascending=False)
    df_grouped["total_ortholog_present"]=df_grouped["module_id"].map(mod_dict)
    return df_grouped
                      


# In[ ]:


#Eukaryotic predicted proteins vs Prokaryotic KEGG module genes


# In[4]:


eukbook_euk_df=output_to_dataframe("eukbookVsKegg.m8")


# In[ ]:


##Filtering and grouping in sample level


# In[41]:


euk_sample_filtered=remove_duplicates(eukbook_euk_df, "sample_id")


# In[74]:


euk_sample_grouped=df_grouper(euk_sample_filtered, "sample_id")


# In[ ]:


##Filtering and grouping in contig level


# In[42]:


euk_contig_filtered=remove_duplicates(eukbook_euk_df, "contig_id")


# In[76]:


euk_contig_grouped=df_grouper(euk_contig_filtered, "contig_id")


# In[ ]:


#Eukaryotic predicted proteins vs Prokaryotic KEGG module genes


# In[40]:


eukbook_unk_df=output_to_dataframe("eukbookVsKegg_unk.m8")


# In[ ]:


##Filtering and grouping in sample level


# In[43]:


unk_sample_filtered=remove_duplicates(eukbook_unk_df, "sample_id")


# In[77]:


unk_sample_grouped=df_grouper(unk_sample_filtered, "sample_id")


# In[ ]:


##Filtering and grouping in contig level


# In[46]:


unk_contig_filtered=remove_duplicates(eukbook_unk_df, "contig_id")


# In[78]:


unk_contig_grouped=df_grouper(unk_contig_filtered, "contig_id")

