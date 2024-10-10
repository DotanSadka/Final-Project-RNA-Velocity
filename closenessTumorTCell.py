import pandas as pd
import numpy as np

def filter_TCell_closeness(is_close_file, file):
    """
    Filter from binary closeness indication in "is_close_file" if a cell was in interaction with tumor
    Save it in a new column in "file"
    """
    is_close = pd.read_csv(is_close_file)
    TCell = pd.read_csv(file)

    filtered_values = is_close[is_close["cell_name"].isin(TCell["Unnamed: 0"])]["count_tumor_total"]

    TCell["isClose"] = list(filtered_values)

    TCell.to_csv(f"df/T-cell_CD8A_{T}.csv", index=False)


if __name__ == "__main__":
    for t in np.arange(0.2, 5.2, 0.2):
        T = round(t,1)
        filter_TCell_closeness("adata_genes_4530.csv", f"../T-cell_CD8A_files/T-cell_CD8A_{T}.csv")
