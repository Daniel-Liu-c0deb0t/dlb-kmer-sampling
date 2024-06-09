import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("distances_closed.csv")
print(df)

df["count"] /= df.groupby(["algorithm", "k"])["count"].transform(sum)
print(df)
ks = sorted(df["k"].unique())

plt.gcf().set_size_inches(7, 5)
plt.suptitle("w = 11", fontsize = 14)

closed_df = df[df["algorithm"] == "Closed syncmer (Edgar 2021)"]
ax = plt.subplot(211)
for k in ks:
    density = closed_df[closed_df["k"] == k]["density"].iloc[0]
    plt.plot(closed_df[closed_df["k"] == k]["distance"], closed_df[closed_df["k"] == k]["count"], label = f"{k} ({density:.3f})")
plt.legend(title = "k (density)", bbox_to_anchor = (1.01, 1))
plt.xlabel("Distance")
plt.ylabel("Relative frequency")
plt.title("Closed syncmer (Edgar 2021)")

optimal_df = df[df["algorithm"] == "Optimal closed syncmer (ours)"]
ax = plt.subplot(212)
for k in ks:
    density = optimal_df[optimal_df["k"] == k]["density"].iloc[0]
    plt.plot(optimal_df[optimal_df["k"] == k]["distance"], optimal_df[optimal_df["k"] == k]["count"], label = f"{k} ({density:.3f})")
plt.legend(title = "k (density)", bbox_to_anchor = (1.01, 1))
plt.xlabel("Distance")
plt.ylabel("Relative frequency")
plt.title("Optimal closed syncmer (ours)")

plt.tight_layout()
plt.savefig("distances_closed.png", bbox_inches = "tight", dpi = 200)
