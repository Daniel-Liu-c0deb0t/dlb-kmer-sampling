import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("metrics.csv")
print(df)

ks = sorted(df["k"].unique())
metrics = df["metric"].unique()
algos = df["algorithm"].unique()

plt.gcf().set_size_inches(9, 5)
plt.suptitle("Density = 1/11", fontsize = 14)

idx = 1
ax_x = None

for i, metric in enumerate(metrics):
    metric_df = df[df["metric"] == metric]
    ax_y = None

    for j, k in enumerate(ks):
        k_df = metric_df[metric_df["k"] == k]
        ax = plt.subplot(len(metrics), len(ks), idx, sharex = ax_x, sharey = ax_y)
        ax_x = ax
        ax_y = ax

        for algo in algos:
            algo_df = k_df[k_df["algorithm"] == algo]
            plt.plot(algo_df["mut_rate"], algo_df["value"] / k_df[k_df["algorithm"] == algos[0]]["value"].values, label = f"{algo}")
            if i == len(metrics) - 1:
                plt.xlabel("Mutation rate")
            if j == 0:
                plt.ylabel(f"Relative {metric}")
            plt.title(f"k = {k}")
        idx += 1

plt.figlegend(algos, title = "Algorithm", loc = "lower right", ncol = len(algos), bbox_to_anchor = (1, -0.1), bbox_transform = plt.gcf().transFigure)
plt.tight_layout()
plt.savefig("metrics.png", bbox_inches = "tight", dpi = 200)
