```
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
%matplotlib inline
```
#read in data
```
motif_position_df = pd.read_csv('/home/jtao/for_shengnan/motif_start_frame_C57BL6J.tsv', sep='\t')
motif_position_df.index = motif_position_df['ID'].values
del motif_position_df['ID']
```
#vehicle position
```
veh_df=motif_position_df[motif_position_df['Factors'].str.contains('c57bl6_atac_veh')]
del veh_df['Factors']
del veh_df['chr']
```
#KLA position 
```
kla_df=motif_position_df[motif_position_df['Factors'].str.contains('c57bl6_atac_kla')]
del kla_df['Factors']
del kla_df['chr']
```
#function to count co-occurence
```
def co_occur(df):
    '''
    input is a dataframe contains all data of the genomic position of motif in chr1.
    output is a dataframe contains the frequency of co-occurance of every two motifs. 
    '''
    motifs = df.columns.values
    count_frame=np.zeros((df.shape[1],df.shape[1]),dtype=np.int)#a dataframe to store future data
    count_frame=pd.DataFrame(count_frame, columns=motifs)
    count_frame.index=motifs
    col_vector=np.zeros((df.shape[0],1),dtype=np.int) #a column vector use in loop
    for i in range (df.shape[1]):#use the column vector to store each column one by one
        col_vector=df.ix[:,i]!= -1# Find if the motifs exist
        for j in range (df.shape[1]): 
            log_v=df.ix[:,j]!= -1#Find if the motifs exist
            vector_sum = 1*col_vector + 1*log_v
            log_input = vector_sum == 2
            count_input=np.sum(1*log_input)#convert logical ou
            count_frame.ix[i,j]=count_frame.ix[i,j]+count_input
            np.fill_diagonal(count_frame.values, 0)# change diagonal back to zero
    return count_frame
```
#count co-occurence
```
veh_cooccurence_df=co_occur(veh_df)
kla_cooccurence_df=co_occur(kla_df)
```
#function to calculate z score of co-occurence
```
def Find_Z(n):
    '''
    For all pairs of motifs - is there a pair that co-occurs more or less often than you would expect.
    input: a dataframe contains frequency of co-occurence of every pair of motifs.
    output:a dataframe contains z score of each pair of motifs.
    '''
    motifs = n.columns.values
    #convert dataframe to matirx
    n=n.as_matrix(columns=None)
    z_matrix = np.zeros((n.shape[0],n.shape[1]-1),dtype=np.float)
    for i in range (n.shape[0]):
        co_motif = n[i,:]
        co_motif = np.delete(co_motif,i)#remove data of the motif co-occur with itself
        z_score=stats.zscore(co_motif)# find z socre 
        z_matrix[i,:]=z_score
    #convert z score matirx to dataframe
    zscore_all=np.zeros((z_matrix.shape[0],z_matrix.shape[0]),dtype=np.float)
    for i in range (z_matrix.shape[0]):
        z_motif_self=z_matrix[i,:]
        z_motif_self=np.insert(z_motif_self,i,100)
        zscore_all[i,:]=z_motif_self
    zscore_frame = pd.DataFrame(zscore_all, columns=motifs)
    zscore_frame.index = motifs
    return zscore_frame
```
#find z scores
```
veh_z=Find_Z(veh_cooccurence_df)
kla_z=Find_Z(kla_cooccurence_df)
```
#function to find co-occurence pdf
```
def Find_pdf(n):
    '''
    For all pairs of motifs - is there a pair that co-occurs more or less often than you would expect.
    input: a matirc contains frequency of co-occurence of every pair of motifs.
    output:an array contains z score of each pair of motifs.
    '''
    motifs = n.columns.values
    #convert dataframe to matirx
    n=n.as_matrix(columns=None)
    p_matrix = np.zeros((n.shape[0],n.shape[1]-1),dtype=np.float)
    for i in range (n.shape[0]):
        #find meand and stv for each motif
        co_motif = n[i,:]
        co_motif = np.delete(co_motif,i)#remove data of the motif co-occur with itself
        mean=np.mean(co_motif)
        std=np.std(co_motif)
        pdf=stats.norm.pdf(co_motif,loc=mean, scale=std)# find z socre 
        p_matrix[i,:]=pdf
    #assign co-occurance z socre of one motif with itself as 100 to make dataframe 
    p_for_dataframe=np.zeros((p_matrix.shape[0],p_matrix.shape[0]),dtype=np.float)
    for i in range (p_matrix.shape[0]):
        p_motif_self=p_matrix[i,:]
        p_motif_self=np.insert(p_motif_self,i,100)
        p_for_dataframe[i,:]=p_motif_self
    p_frame = pd.DataFrame(p_for_dataframe, columns=motifs)
    p_frame.index = motifs
    return p_frame
```
#find pdf
```
veh_cooccurence_pdf=Find_pdf(veh_cooccurence_df)*195
kla_cooccurence_pdf=Find_pdf(kla_cooccurence_df)*195
```
