from dask.distributed import Client, Future, as_completed
import dask.dataframe as dd
import pandas as pd
from pandas.errors import EmptyDataError

def process_file(file_path: str, snps_to_filter: Future) -> dd.DataFrame:
    try:
        ddf = dd.read_csv(
            file_path.strip(),
            usecols=['Chr', 'SNP', 'bp', 'A1', 'A2', 'b', 'se'],
            dtype={'SNP': 'object'},
            sep="\t"
        )
        ddf = ddf[ddf['SNP'].isin(snps_to_filter)]
        ddf['z'] = ddf['b'] / ddf['se']
        return ddf
    except EmptyDataError:
        return None

def main():
    # Start Dask Client
    client = Client(n_workers=8, threads_per_worker=1)

    # Read SNP list and scatter it across the cluster
    snps_to_filter = pd.read_csv("/seq/vgb/dd/gwas/geno/20230905.prune.in", header=None)[0].values.tolist()
    snps_future = client.scatter(snps_to_filter, broadcast=True)
    
    # Read file list once and store it
    file_list = open("/seq/vgb/dd/gwas/assoc/gwa-to-mat_thesis.txt").readlines()
    
    # Initialize list to hold Dask Futures
    futures_list = []
    
    for i, file_path in enumerate(file_list):
        future = client.submit(process_file, file_path, snps_future)
        futures_list.append(future)

    # Gather results as they come in
    ddf_list = []
    for future in as_completed(futures_list):
        result = future.result()
        if result is not None:
            result['phenotype'] = i
            ddf_list.append(result)
    
    # Concatenate all Dask DataFrames
    final_ddf = dd.concat(ddf_list)
    
    # Compute and convert to Pandas DataFrame
    final_df = final_ddf.compute()
    
    # Rest of the code remains identical for saving CSVs and pivoting

if __name__ == '__main__':
    main()