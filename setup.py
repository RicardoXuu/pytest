from setuptools import setup

setup(
    name='SingleCellGGM',
    version='0.1',
    license='BSD-3-Clause',
    author='MaShisongLab',
    py_modules=['SingleCellGGM'],
    description='An algorithm to conduct single-cell gene co-expression network analysis',
    long_description='SingleCellGGM is an algorithm to conduct gene co-expression network analysis for single-cell transcriptome datasets based on the graphical Gaussian model (GGM). It uses a process consisting of multiple iterations to calculate partial correlation coefficients (pcors) between genes for co-expression network construction. For more details about the algorithm, please refer to Xu et. al 2023 (https://doi.org/10.1101/2023.02.05.526424). This module contains a class, also called SingleCellGGM, for conducting co-expression network analysis. It takes a log-normalized single-cell gene expression matrix (samples in rows and genes in columns), the number of iterations (e.g. 2000 or 20000), an array containing the gene names, and the name of dataset (optional) as inputs. The format is: ggm = SingleCellGGM(expression_matrix, number_of_iterations,  gene_names,  dataset_name). There is another function, called getsigedges, under this class. After building object named ggm based on input data, one can filter significant edges according to different thresholds by running the getsigedges function, and the results will be updated to ggm.SigEdges. The format is: ggm.getsigedges(pcor_threshold, coex_cell_threshold), where pcor_threshold represents minimum threshold of pcor of reserved gene pair; coex_cell_threshold represents minimum threshold of co-expressed cells of reserved gene pair.',
    url='https://github.com/RicardoXuu/pytest',
    keywords='SingleCellGGM',
    install_requires=['numpy', 'pandas', 'numba'],
)
