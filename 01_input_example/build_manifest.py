import pandas as pd

df = pd.read_csv("manifest.tsv", sep="\t")

fwd = df[df["direction"] == "forward"][["sample-id", "absolute-filepath"]]
fwd.columns = ["sample-id", "forward-absolute-filepath"]

rev = df[df["direction"] == "reverse"][["sample-id", "absolute-filepath"]]
rev.columns = ["sample-id", "reverse-absolute-filepath"]

manifest_v2 = pd.merge(fwd, rev, on="sample-id")
manifest_v2.to_csv("manifest_v2.tsv", sep="\t", index=False)

print("manifest_v2.tsv created successfully")
