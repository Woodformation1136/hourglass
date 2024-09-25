import multiprocessing as mp
import pandas as pd
import numpy as np
import sys

# get sys args
argv = sys.argv[1:]
if len(argv) != 2:
    print("Usage: python create_pairwise_distance_matrix.py <projection.csv> <output.csv>")
    sys.exit(1)

projection_csv, output_csv = argv

# read projection.csv
projections = pd.read_csv(projection_csv, index_col=0)
n_rows, n_cols = projections.shape
    
# create pairwise distance matrix
pairwise_distance_matrix = pd.DataFrame(index=projections.index, columns=projections.index)

# report start time
print("Calculating pairwise distance matrix...")
print("Number of rows: ", n_rows)
print("Start time: ", pd.Timestamp.now())
for i in range(n_rows):
    for j in range(i, n_rows):
        distance = np.sqrt(np.sum(np.square(projections.iloc[i] - projections.iloc[j])))
        pairwise_distance_matrix.iloc[i, j] = distance
        pairwise_distance_matrix.iloc[j, i] = distance
    print("[%s] Progress: %i/%i" % (pd.Timestamp.now(), i+1, n_rows))
print("Finished time: ", pd.Timestamp.now())

# save pairwise distance matrix
pairwise_distance_matrix.to_csv(output_csv)
