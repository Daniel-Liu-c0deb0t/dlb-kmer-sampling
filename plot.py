import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results.csv")
print(df)

df["count"] /= df.groupby(["algorithm", "k"])["count"].transform(sum)
print(df)
ks = sorted(df["k"].unique())

plt.gcf().set_size_inches(7, 7)
plt.suptitle("Density = 1/11", fontsize = 14)

start_df = df[df["algorithm"] == "Open syncmer (Edgar 2021)"]
ax = plt.subplot(311)
for k in ks:
    plt.plot(start_df[start_df["k"] == k]["distance"], start_df[start_df["k"] == k]["count"], label = f"{k}")
plt.legend(title = "k")
min_dist = min(start_df["distance"])
plt.axvline(min_dist, ls = ":", lw = 1, c = "C0")
plt.text(min_dist, -0.02, f"{min_dist}", ha = "center", va = "top", fontsize = 8, c = "C0", transform = plt.gca().get_xaxis_transform())
plt.xlabel("Distance")
plt.ylabel("Relative frequency")
plt.title("Open syncmer (Edgar 2021)")

mid_df = df[df["algorithm"] == "Open syncmer (Shaw & Yu 2022)"]
plt.subplot(312, sharex = ax)
for k in ks:
    plt.plot(mid_df[mid_df["k"] == k]["distance"], mid_df[mid_df["k"] == k]["count"], label = f"{k}")
plt.legend(title = "k")
min_dist = min(mid_df["distance"])
plt.axvline(min_dist, ls = ":", lw = 1, c = "C0")
plt.text(min_dist, -0.02, f"{min_dist}", ha = "center", va = "top", fontsize = 8, c = "C0", transform = plt.gca().get_xaxis_transform())
plt.xlabel("Distance")
plt.ylabel("Relative frequency")
plt.title("Open syncmer (Shaw & Yu 2022)")

dlb_df = df[df["algorithm"] == "DLB (ours)"]
plt.subplot(313, sharex = ax)
for i, k in enumerate(ks):
    plt.plot(dlb_df[dlb_df["k"] == k]["distance"], dlb_df[dlb_df["k"] == k]["count"], label = f"{k}")
    min_dist = min(dlb_df[dlb_df["k"] == k]["distance"])
    plt.axvline(min_dist, ls = ":", lw = 1, c = f"C{i}")
    plt.text(min_dist, -0.02, f"{min_dist}", ha = "center", va = "top", fontsize = 8, c = f"C{i}", transform = plt.gca().get_xaxis_transform())
plt.legend(title = "k")
plt.xlabel("Distance")
plt.ylabel("Relative frequency")
plt.title("DLB (ours)")

plt.tight_layout()
plt.savefig("distances.png", bbox_inches = "tight", dpi = 200)
