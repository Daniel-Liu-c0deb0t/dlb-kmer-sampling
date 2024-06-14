import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("density_mod.csv")
print(df)

algos = df["algorithm"].unique()

w = 31

plt.gcf().set_size_inches(7, 5)
plt.suptitle(f"w = {w}", fontsize = 14)

ax = plt.subplot(211)
for algo in algos:
    algo_df = df[df["algorithm"] == algo]
    plt.plot(algo_df["k"], algo_df["density"], label = f"{algo}")
plt.axvline(w, ls = ":", lw = 1, c = "red")
plt.text(w + 0.3, 0.2, f"k = {w}", ha = "left", va = "bottom", c = "red", transform = plt.gca().get_xaxis_transform())
rand = 2 / (w + 1)
plt.axhline(rand, lw = 1, c = "black")
plt.text(0.005, rand - 0.0003, "Random", ha = "left", va = "top", c = "black", transform = plt.gca().get_yaxis_transform())
lb = 1 / w
plt.axhline(lb, lw = 1, c = "black")
plt.text(0.005, lb, "Lower bound", ha = "left", va = "bottom", c = "black", transform = plt.gca().get_yaxis_transform())
plt.legend(title = "Algorithm")
plt.xlabel("k")
plt.ylabel("Density")

ax = plt.subplot(212)
for algo in algos:
    algo_df = df[df["algorithm"] == algo]
    plt.plot(algo_df["k"], algo_df["w"], label = f"{algo}")
plt.legend(title = "Algorithm")
plt.xlabel("k")
plt.ylabel("w' (w' <= w)")

plt.tight_layout()
plt.savefig("density_mod.png", bbox_inches = "tight", dpi = 200)
