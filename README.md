```python
# Import python modules
import scanpy as sc
import numpy as np
import scipy.sparse as sp
```


```python
# Read in the single-cell gene expression matrix as adata by scanpy
adata = sc.read_10x_mtx('./pbmc3k/filtered_gene_bc_matrices/hg19',var_names='gene_ids',cache=True)
sc.pp.filter_genes(adata,min_cells=10)
```


```python
adata.X[1:5,1:5]
```




    <4x4 sparse matrix of type '<class 'numpy.float32'>'
    	with 2 stored elements in Compressed Sparse Row format>




```python
# Type1: Get gene expression matrix (dense matrix)
#expression_matrix = adata.X.todense()
#expression_matrix = np.transpose(expression_matrix)
# Ensure the matrix is in double format
#expression_matrix = np.float64(expression_matrix)
```


```python
# Type2: Get gene expression matrix (ndarray matrix)
expression_matrix = adata.X.toarray()
expression_matrix = np.transpose(expression_matrix)
# Ensure the matrix is in double format
expression_matrix = np.float64(expression_matrix)
```


```python
# Type3: Get gene expression matrix (sparse csr matrix)
#expression_matrix = adata.X
#expression_matrix = np.transpose(expression_matrix)
# Ensure the matrix is in double format
#expression_matrix = expression_matrix.astype('float64')
```


```python
print(type(expression_matrix))
expression_matrix.shape
```

    <class 'numpy.ndarray'>





    (11139, 2700)




```python
# Log normalization
expression_matrix = np.log2((expression_matrix / np.sum(expression_matrix, axis=0)) * 10000 + 1)
```


```python
expression_matrix = np.transpose(expression_matrix)
```


```python
print(type(expression_matrix))
expression_matrix.shape
```

    <class 'numpy.ndarray'>





    (2700, 11139)




```python
expression_matrix[0:5,0:5]
```




    array([[0.        , 0.        , 0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ],
           [0.        , 0.        , 0.        , 2.06373126, 0.        ],
           [0.        , 0.        , 0.        , 5.13727927, 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ]])




```python
# Get gene names as an array
gene = adata.var.index
gene = np.array(gene)
```


```python
# Import SingleCellGGM modules
from SingleCellGGM import *
```


```python
help(SingleCellGGM)
```

    Help on class SingleCellGGM in module SingleCellGGM:
    
    class SingleCellGGM(builtins.object)
     |  SingleCellGGM(x, round_num, gene_name, dataset_name='na')
     |  
     |  Methods defined here:
     |  
     |  __init__(self, x, round_num, gene_name, dataset_name='na')
     |      x : an expression matrix, cells in row, gene in column; Test3
     |      round_num : number of iterations;
     |      gene_name : an array of gene names; 
     |      dataset_name : optional
     |  
     |  adjust_cutoff(self, pcor_threshold=0.03, coex_cell_threshold=10)
     |      pcor_threshold : minimum threshold of pcor of reserved gene pair;
     |      coex_cell_threshold : minimum threshold of co-expressed cells of reserved gene pair
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  RoundNumber = []
     |  
     |  SigEdges = []
     |  
     |  coexpressed_cell_num = []
     |  
     |  data_name = []
     |  
     |  gene_name = []
     |  
     |  gene_num = []
     |  
     |  pcor_all = []
     |  
     |  pcor_sampling_num = []
     |  
     |  rho_all = []
     |  
     |  samples_num = []
    



```python
help(SingleCellGGM.adjust_cutoff)
```

    Help on function adjust_cutoff in module SingleCellGGM:
    
    adjust_cutoff(self, pcor_threshold=0.03, coex_cell_threshold=10)
        pcor_threshold : minimum threshold of pcor of reserved gene pair;
        coex_cell_threshold : minimum threshold of co-expressed cells of reserved gene pair
    



```python
help(fdr_control)
```

    Help on class fdr_control in module SingleCellGGM:
    
    class fdr_control(builtins.object)
     |  fdr_control(x, ggm_ori, permutation_fraction=1)
     |  
     |  Methods defined here:
     |  
     |  __init__(self, x, ggm_ori, permutation_fraction=1)
     |      x : gene expression matrix; Test6
     |      ggm_ori : result from SingleCellGGM function;
     |      permutation_fraction : fraction of genes to be permutated (default is 1);
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  RoundNumber = []
     |  
     |  SigEdges = []
     |  
     |  data_name = []
     |  
     |  fdr = []
     |  
     |  gene_name = []
     |  
     |  gene_num = []
     |  
     |  samples_num = []
    



```python
# Conduct single-cell gene co-expression analysis via SingleCellGGM
ggm = SingleCellGGM( expression_matrix, 20000, gene, "pbmc3k" )
```

    Calculating pcor in 20000 iterations.
    
    Estimated to complete in 6.86 hours.
    
    1000 iterations done. 1.01 sec/iteration in average. 5.35 hours to go.
    
    2000 iterations done. 1.02 sec/iteration in average. 5.09 hours to go.
    
    3000 iterations done. 1.03 sec/iteration in average. 4.84 hours to go.
    
    4000 iterations done. 1.10 sec/iteration in average. 4.87 hours to go.
    
    5000 iterations done. 0.98 sec/iteration in average. 4.07 hours to go.
    
    6000 iterations done. 1.06 sec/iteration in average. 4.12 hours to go.
    
    7000 iterations done. 0.99 sec/iteration in average. 3.57 hours to go.
    
    8000 iterations done. 0.99 sec/iteration in average. 3.30 hours to go.
    
    9000 iterations done. 0.99 sec/iteration in average. 3.03 hours to go.
    
    10000 iterations done. 1.01 sec/iteration in average. 2.82 hours to go.
    
    11000 iterations done. 0.99 sec/iteration in average. 2.49 hours to go.
    
    12000 iterations done. 1.00 sec/iteration in average. 2.21 hours to go.
    
    13000 iterations done. 1.00 sec/iteration in average. 1.94 hours to go.
    
    14000 iterations done. 1.00 sec/iteration in average. 1.67 hours to go.
    
    15000 iterations done. 1.00 sec/iteration in average. 1.39 hours to go.
    
    16000 iterations done. 1.00 sec/iteration in average. 1.12 hours to go.
    
    17000 iterations done. 1.01 sec/iteration in average. 0.84 hours to go.
    
    18000 iterations done. 1.01 sec/iteration in average. 0.56 hours to go.
    
    19000 iterations done. 1.01 sec/iteration in average. 0.28 hours to go.
    
    20000 iterations done. 1.01 sec/iteration in average. 0.00 hours to go.
    
    20000 iterations done.
    



```python
# Examine the results
ggm.SigEdges
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GeneA</th>
      <th>GeneB</th>
      <th>Pcor</th>
      <th>SamplingTime</th>
      <th>r</th>
      <th>Cell_num_A</th>
      <th>Cell_num_B</th>
      <th>Cell_num_coexpressed</th>
      <th>Dataset</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000008128</td>
      <td>ENSG00000248333</td>
      <td>0.461646</td>
      <td>641</td>
      <td>0.530299</td>
      <td>164</td>
      <td>62</td>
      <td>50</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000173369</td>
      <td>ENSG00000173372</td>
      <td>0.674732</td>
      <td>656</td>
      <td>0.736493</td>
      <td>17</td>
      <td>36</td>
      <td>14</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000126709</td>
      <td>ENSG00000187608</td>
      <td>0.036944</td>
      <td>620</td>
      <td>0.372184</td>
      <td>1082</td>
      <td>1206</td>
      <td>671</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000223759</td>
      <td>ENSG00000235999</td>
      <td>0.572937</td>
      <td>644</td>
      <td>0.645156</td>
      <td>31</td>
      <td>56</td>
      <td>26</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000150337</td>
      <td>ENSG00000198019</td>
      <td>0.122937</td>
      <td>663</td>
      <td>0.289769</td>
      <td>119</td>
      <td>61</td>
      <td>25</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>397</th>
      <td>ENSG00000198727</td>
      <td>ENSG00000198763</td>
      <td>0.088896</td>
      <td>606</td>
      <td>0.312927</td>
      <td>2517</td>
      <td>2416</td>
      <td>2293</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>398</th>
      <td>ENSG00000198727</td>
      <td>ENSG00000198804</td>
      <td>0.079614</td>
      <td>606</td>
      <td>0.202792</td>
      <td>2517</td>
      <td>2686</td>
      <td>2511</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>399</th>
      <td>ENSG00000198727</td>
      <td>ENSG00000198899</td>
      <td>0.069512</td>
      <td>625</td>
      <td>0.265027</td>
      <td>2517</td>
      <td>2014</td>
      <td>1925</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>400</th>
      <td>ENSG00000198727</td>
      <td>ENSG00000198886</td>
      <td>0.052989</td>
      <td>629</td>
      <td>0.280573</td>
      <td>2517</td>
      <td>2588</td>
      <td>2433</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>401</th>
      <td>ENSG00000220023</td>
      <td>ENSG00000149531</td>
      <td>0.427986</td>
      <td>614</td>
      <td>0.506088</td>
      <td>323</td>
      <td>277</td>
      <td>165</td>
      <td>pbmc3k</td>
    </tr>
  </tbody>
</table>
<p>402 rows × 9 columns</p>
</div>




```python
# Save all gene pairs to a file for gene co-expression network construction
ggm.SigEdges.iloc[:, :3].to_csv('pbmc3k.ggm.network.txt', sep='\t', index=False, header=True)
```


```python
# FDR Control - control the false discovery rate with different pcor cutoff value
fdr = fdr_control(expression_matrix, ggm)
```

    Calculating pcor in 20000 iterations.
    
    Estimated to complete in 6.18 hours.
    
    1000 iterations done. 1.01 sec/iteration in average. 5.34 hours to go.
    
    2000 iterations done. 1.00 sec/iteration in average. 5.02 hours to go.
    
    3000 iterations done. 1.01 sec/iteration in average. 4.76 hours to go.
    
    4000 iterations done. 1.01 sec/iteration in average. 4.49 hours to go.
    
    5000 iterations done. 1.01 sec/iteration in average. 4.23 hours to go.
    
    6000 iterations done. 1.02 sec/iteration in average. 3.97 hours to go.
    
    7000 iterations done. 1.03 sec/iteration in average. 3.71 hours to go.
    
    8000 iterations done. 1.02 sec/iteration in average. 3.40 hours to go.
    
    9000 iterations done. 1.03 sec/iteration in average. 3.16 hours to go.
    
    10000 iterations done. 1.03 sec/iteration in average. 2.85 hours to go.
    
    11000 iterations done. 1.03 sec/iteration in average. 2.59 hours to go.
    
    12000 iterations done. 1.03 sec/iteration in average. 2.30 hours to go.
    
    13000 iterations done. 1.03 sec/iteration in average. 2.01 hours to go.
    
    14000 iterations done. 1.04 sec/iteration in average. 1.73 hours to go.
    
    15000 iterations done. 1.04 sec/iteration in average. 1.45 hours to go.
    
    16000 iterations done. 1.20 sec/iteration in average. 1.33 hours to go.
    
    17000 iterations done. 1.21 sec/iteration in average. 1.01 hours to go.
    
    18000 iterations done. 1.48 sec/iteration in average. 0.82 hours to go.
    
    19000 iterations done. 1.18 sec/iteration in average. 0.33 hours to go.
    
    20000 iterations done. 1.13 sec/iteration in average. 0.00 hours to go.
    
    20000 iterations done.
    



```python
fdr.fdr
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Pcor</th>
      <th>SigEdgeNum</th>
      <th>FDR</th>
      <th>SigPermutatedEdgeNum</th>
      <th>SigPermutatedEdgeProportion</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.010</td>
      <td>1463</td>
      <td>0.327</td>
      <td>478</td>
      <td>7.71e-06</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.011</td>
      <td>1318</td>
      <td>0.313</td>
      <td>413</td>
      <td>6.66e-06</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.012</td>
      <td>1212</td>
      <td>0.291</td>
      <td>353</td>
      <td>5.69e-06</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.013</td>
      <td>1109</td>
      <td>0.264</td>
      <td>293</td>
      <td>4.72e-06</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.014</td>
      <td>1019</td>
      <td>0.242</td>
      <td>247</td>
      <td>3.98e-06</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>86</th>
      <td>0.096</td>
      <td>87</td>
      <td>&lt; 0.00662</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>87</th>
      <td>0.097</td>
      <td>86</td>
      <td>&lt; 0.00662</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>88</th>
      <td>0.098</td>
      <td>86</td>
      <td>&lt; 0.00662</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>89</th>
      <td>0.099</td>
      <td>84</td>
      <td>&lt; 0.00662</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>90</th>
      <td>0.100</td>
      <td>83</td>
      <td>&lt; 0.00662</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>91 rows × 5 columns</p>
</div>




```python
# Adjust pcor cutoff and coexpressed cell cutoff
ggm.adjust_cutoff(0.02, 10)
ggm.SigEdges
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GeneA</th>
      <th>GeneB</th>
      <th>Pcor</th>
      <th>SamplingTime</th>
      <th>r</th>
      <th>Cell_num_A</th>
      <th>Cell_num_B</th>
      <th>Cell_num_coexpressed</th>
      <th>Dataset</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>P_ENSG00000186827</td>
      <td>P_ENSG00000186891</td>
      <td>0.023259</td>
      <td>620</td>
      <td>0.194450</td>
      <td>155</td>
      <td>92</td>
      <td>29</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>1</th>
      <td>P_ENSG00000008128</td>
      <td>P_ENSG00000248333</td>
      <td>0.461646</td>
      <td>641</td>
      <td>0.530299</td>
      <td>164</td>
      <td>62</td>
      <td>50</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>2</th>
      <td>P_ENSG00000173369</td>
      <td>P_ENSG00000173372</td>
      <td>0.674732</td>
      <td>656</td>
      <td>0.736493</td>
      <td>17</td>
      <td>36</td>
      <td>14</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>3</th>
      <td>P_ENSG00000011009</td>
      <td>P_ENSG00000049239</td>
      <td>0.024675</td>
      <td>634</td>
      <td>0.111930</td>
      <td>231</td>
      <td>51</td>
      <td>14</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>4</th>
      <td>P_ENSG00000126709</td>
      <td>P_ENSG00000187608</td>
      <td>0.036944</td>
      <td>620</td>
      <td>0.372184</td>
      <td>1082</td>
      <td>1206</td>
      <td>671</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>707</th>
      <td>P_ENSG00000198727</td>
      <td>P_ENSG00000198804</td>
      <td>0.079614</td>
      <td>606</td>
      <td>0.202792</td>
      <td>2517</td>
      <td>2686</td>
      <td>2511</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>708</th>
      <td>P_ENSG00000198727</td>
      <td>P_ENSG00000198712</td>
      <td>0.029128</td>
      <td>677</td>
      <td>0.269995</td>
      <td>2517</td>
      <td>2460</td>
      <td>2318</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>709</th>
      <td>P_ENSG00000198727</td>
      <td>P_ENSG00000198899</td>
      <td>0.069512</td>
      <td>625</td>
      <td>0.265027</td>
      <td>2517</td>
      <td>2014</td>
      <td>1925</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>710</th>
      <td>P_ENSG00000198727</td>
      <td>P_ENSG00000198886</td>
      <td>0.052989</td>
      <td>629</td>
      <td>0.280573</td>
      <td>2517</td>
      <td>2588</td>
      <td>2433</td>
      <td>pbmc3k</td>
    </tr>
    <tr>
      <th>711</th>
      <td>P_ENSG00000220023</td>
      <td>P_ENSG00000149531</td>
      <td>0.427986</td>
      <td>614</td>
      <td>0.506088</td>
      <td>323</td>
      <td>277</td>
      <td>165</td>
      <td>pbmc3k</td>
    </tr>
  </tbody>
</table>
<p>712 rows × 9 columns</p>
</div>




```python
fdr.fdr.to_csv('pbmc3k.ggm.fdr.txt',sep='\t', index = False, header = True)
```


```python

```
