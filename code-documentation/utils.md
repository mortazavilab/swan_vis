# Module swan\_vis.utils

Functions
---------


`calc_pi(adata, t_df, obs_col='dataset')`
:   Calculate the percent isoform per gene per condition given by `obs_col`.
    Default column to use is `adata.obs` index column, `dataset`.

    Parameters:
            adata (anndata AnnData): Annotated data object from the SwanGraph
            t_df (pandas DataFrame): Pandas Dataframe that has index to
                    gene id mapping
            obs_col (str): Column name from adata.obs table to group on.
                    Default: 'dataset'

    Returns:
            df (pandas DataFrame): Pandas DataFrame where rows are the different
                    conditions from `obs_col` and the columns are transcript ids in the
                    SwanGraph, and values represent the percent isoform usage per gene
                    per condition.
            sums (pandas DataFrame): Pandas DataFrame where rows are the different
                    conditions from `obs_col` and the columns are transcript ids in the
                    SwanGraph, and values represent the cumulative counts per isoform
                    per condition.


`calc_total_counts(adata, obs_col='dataset', layer='counts')`
:   Calculate cumulative expression per adata entry based on condition given
    by `obs_col`. Default column to use is `adata.obs` index column, `dataset`.

    Parameters:
            adata (anndata AnnData): Annotated data object from the SwanGraph
            obs_col (str): Column name from adata.obs table to group on.
                    Default: 'dataset'
            layer (str): Layer of AnnData to pull from. Default = 'counts'

    Returns:
            df (pandas DataFrame): Pandas DataFrame where rows are the different
                    conditions from `obs_col` and the columns are transcript ids in the
                    SwanGraph, and values represent the cumulative counts per isoform
                    per condition.


`calc_tpm(adata, sg_df=None, obs_col='dataset')`
:   Calculate the TPM per condition given by `obs_col`.
    Default column to use is `adata.obs` index column, `dataset`.

    Parameters:
            adata (anndata AnnData): Annotated data object from the SwanGraph
            sg_df (pandas DataFrame): Pandas DataFrame from SwanGraph that will
                    be used to order the rows of resultant TPM DataFrame
            obs_col (str or list of str): Column name from adata.obs table to group on.
                    Default: 'dataset'

    Returns:
            df (pandas DataFrame): Pandas datafrom where rows are the different
                    conditions from `obs_col` and the columns are transcript ids in the
                    SwanGraph, and values represent the TPM value per isoform per
                    condition.

`read(file)`
:   Read a SwanGraph from a saved pickle file.

    Parameters:
      file (str): Name / path to saved Swan file

    Returns:
      sg (SwanGraph): SwanGraph stored in file


`save_fig(oname)`
:   Save the current figure as a png with a given file prefix.

  	Parameters:
  		oname (str): Path / prefix to saved image


`validate_gtf(fname)`
:	Validates that the input GTF is valid input to Swan.

	Parameters:
		fname (str): Path / name of GTF file
