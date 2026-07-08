import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova

PANINI_FILE = "panini-pirate2Rtab.Rtab.csv"
META_FILE = "microreact-2.csv"

panini = pd.read_csv(PANINI_FILE)
meta = pd.read_csv(META_FILE)

# PANINI id looks like 105577_PA010; metadata id is 105577
panini["merge_id"] = panini["id"].astype(str).str.extract(r"^(\d+)")[0].astype(int)
meta["merge_id"] = meta["id"].astype(int)

df = panini.merge(meta, on="merge_id", how="inner", suffixes=("_panini", "_meta"))

df = df.dropna(subset=["x", "y", "source_label", "clonal_complex (MLST)"])

print(f"PANINI rows: {len(panini)}")
print(f"Metadata rows: {len(meta)}")
print(f"Merged rows used: {len(df)}")

coords = df[["x", "y"]].to_numpy()

# Distance matrix from PANINI coordinates
dist = squareform(pdist(coords, metric="euclidean"))
ids = df["id_panini"].astype(str).tolist()
dm = DistanceMatrix(dist, ids=ids)

def r2_from_coords(coords, groups):
    groups = pd.Series(groups).astype(str).values
    grand = coords.mean(axis=0)
    ss_total = ((coords - grand) ** 2).sum()
    ss_between = 0.0
    for g in np.unique(groups):
        sub = coords[groups == g]
        centroid = sub.mean(axis=0)
        ss_between += len(sub) * ((centroid - grand) ** 2).sum()
    return ss_between / ss_total

tests = {
    "host_source": df["source_label"].astype(str).values,
    "clonal_complex": df["clonal_complex (MLST)"].astype(str).values,
}

with open("PANINI_PERMANOVA_results.txt", "w") as out:
    for name, groups in tests.items():
        print(f"\nPERMANOVA: {name}")
        res = permanova(dm, groups, permutations=999)
        r2 = r2_from_coords(coords, groups)

        print(res)
        print(f"Approximate R2 from PANINI coordinates: {r2:.4f}")

        out.write(f"\nPERMANOVA: {name}\n")
        out.write(str(res))
        out.write(f"\nApproximate R2 from PANINI coordinates: {r2:.4f}\n")

print("\nDone. Results written to PANINI_PERMANOVA_results.txt")