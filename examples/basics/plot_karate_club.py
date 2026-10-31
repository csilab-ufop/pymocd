"""
# Communities in Zachary's karate club

Run [`pymocd.rimpso`][pymocd.rimpso] on Zachary's karate club, draw the
partition it returns, and score it against the club's recorded split.
"""

# %%
# ## Detect communities
#
# Every detector takes a NetworkX graph with integer node ids (or an igraph
# graph) and returns a `dict` mapping each node to a community id.

import matplotlib.pyplot as plt
import networkx as nx

import pymocd

G = nx.karate_club_graph()
partition = pymocd.rimpso(G)

print(f"{G.number_of_nodes()} nodes, {len(set(partition.values()))} communities")

# %%
# ## Draw the partition
#
# One colour per community. Isolated nodes carry community `-1` and would be
# drawn as a group of their own; drop them from the dict first if you do not
# want them shown.

pos = nx.spring_layout(G, seed=7)
communities = sorted(set(partition.values()))
cmap = plt.colormaps["tab10"] if len(communities) <= 10 else plt.colormaps["tab20"].resampled(len(communities))

fig, ax = plt.subplots(figsize=(6, 4.5))
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.25)
for i, c in enumerate(communities):
    nodes = [v for v in G if partition[v] == c]
    nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=[cmap(i)], node_size=180, ax=ax, label=f"community {c}")
nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
ax.legend(loc="lower left", fontsize=8, frameon=False)
ax.set_axis_off()
fig.tight_layout()
plt.show()

# %%
# ## Score against ground truth
#
# NetworkX records which faction each member joined after the split in the
# `club` node attribute. [`pymocd.gt_metrics`][pymocd.gt_metrics] returns NMI,
# AMI, ARI and pairwise F1 in one call.

gt = {v: int(G.nodes[v]["club"] != "Mr. Hi") for v in G}
nmi, ami, ari, f1 = pymocd.gt_metrics(partition, gt)
print(f"NMI={nmi:.3f}  AMI={ami:.3f}  ARI={ari:.3f}  F1={f1:.3f}")
