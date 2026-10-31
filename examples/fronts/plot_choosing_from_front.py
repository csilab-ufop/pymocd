"""
# Choosing a partition from a front

Take HP-MOCD's Pareto front on Zachary's karate club from
[`pymocd.hpmocd_fronts`][pymocd.hpmocd_fronts] and pick one member in four
ways: three selection rules (maximum modularity, a target number of
communities, and RIMPSO's label-free rule through
[`pymocd.rimpso_select`][pymocd.rimpso_select]) and the ground-truth oracle.
"""

# %%
# ## The front
#
# [`pymocd.hpmocd_fronts`][pymocd.hpmocd_fronts] returns HP-MOCD's Pareto
# front, the candidate set that [`pymocd.hpmocd`][pymocd.hpmocd] selects from,
# as a list of partitions. A front can list the same partition more than once,
# so the code keeps one copy of each. HP-MOCD has no seed parameter, and its
# module README records that two calls on the same graph differ, so the front
# and every number on this page can change from one build of the documentation
# to the next.

from collections import Counter

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import pymocd

G = nx.karate_club_graph()
truth = {v: int(G.nodes[v]["club"] != "Mr. Hi") for v in G}  # the recorded split


def as_sets(partition):
    """dict[node, community] -> list of node sets, the form NetworkX expects."""
    groups = {}
    for v, c in partition.items():
        groups.setdefault(c, set()).add(v)
    return list(groups.values())


front, seen = [], set()
for p in pymocd.hpmocd_fronts(G):
    key = frozenset(map(frozenset, as_sets(p)))
    if key not in seen:
        seen.add(key)
        front.append(p)

k = np.array([len(set(p.values())) for p in front])
print(f"{len(front)} distinct partitions on the front, from {k.min()} to {k.max()} communities")

# %%
# ## The objectives of each member
#
# HP-MOCD minimizes the two objectives of
# [Shi et al. (2012)](https://doi.org/10.1016/j.asoc.2011.10.005):
# $\mathrm{intra}$, the fraction of edges that run between communities, and
# $\mathrm{inter} = \sum_c (d_c / 2m)^2$, where $d_c$ is the degree sum of
# community $c$ and $m$ the number of edges. Their sum is $1 - Q$, with $Q$ the
# modularity ([Core concepts](../../guide/concepts.md)). `hpmocd_fronts` returns
# partitions only, so the code computes both objectives from these
# definitions; the [`pymocd.HpMocd`][pymocd.HpMocd] class returns them with
# each member through
# [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front], from a
# separate run. pymocd ignores edge weights, and NetworkX's karate graph has
# them, so modularity is computed with `weight=None`.

m = G.number_of_edges()


def intra_inter(partition):
    internal = sum(partition[u] == partition[v] for u, v in G.edges())
    degree_sum = Counter()
    for v, d in G.degree():
        degree_sum[partition[v]] += d
    return 1 - internal / m, sum((d / (2 * m)) ** 2 for d in degree_sum.values())


intra, inter = np.array([intra_inter(p) for p in front]).T
Q = np.array([nx.community.modularity(G, as_sets(p), weight=None) for p in front])
ami = np.array([pymocd.ami(p, truth) for p in front])
print(f"largest |Q - (1 - intra - inter)| over the front: {np.abs(Q - (1 - intra - inter)).max():.1e}")

# %%
# ## Four ways to choose
#
# - **Maximum modularity.** The member with the highest $Q$, the rule that
#   `hpmocd` applies to its own front.
# - **Closest to a target number of communities.** The member whose number of
#   communities is closest to a number fixed in advance, here 2, the number of
#   factions the club split into; ties go to the higher $Q$. The rule needs that
#   number.
# - **RIMPSO's rule.** [`pymocd.rimpso_select`][pymocd.rimpso_select] applies
#   RIMPSO's label-free selector, a degree-corrected assortative block model
#   with a Schwarz penalty, to any list of candidate partitions. It returns a two-element list:
#   the selected index and the `[cut, pair]` point of every candidate.
# - **Front oracle.** The member with the highest AMI against the recorded
#   split. It reads the ground truth, so it is a yardstick for evaluation, not
#   a rule that can be applied to a graph without labels.

target_k = 2
selected_index, cut_pair = pymocd.rimpso_select(G, front)

choices = {
    "maximum modularity": int(np.argmax(Q)),
    f"closest to k = {target_k}": min(range(len(front)), key=lambda i: (abs(k[i] - target_k), -Q[i])),
    "RIMPSO's rule": selected_index,
    "oracle (uses labels)": int(np.argmax(ami)),
}
print(f"rimpso_select: index {selected_index}, {len(cut_pair)} [cut, pair] points, "
      f"selected point {[round(x, 4) for x in cut_pair[selected_index]]}")

BLUE, ORANGE, AQUA, INK = "#2a78d6", "#eb6834", "#1baf7a", "#1a1a19"
style = {  # rule -> (color, ring size); rings of different sizes stay visible when choices coincide
    "maximum modularity": (BLUE, 420),
    f"closest to k = {target_k}": (ORANGE, 250),
    "RIMPSO's rule": (AQUA, 120),
}

fig, (ax_front, ax_bar) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})

ax_front.scatter(intra, inter, s=22, color="#9a9994", zorder=2, label="front member")
q_best = Q[choices["maximum modularity"]]
line = np.linspace(0, 1, 2)
ax_front.plot(line, 1 - q_best - line, color=BLUE, lw=1, ls="--", zorder=1,
              label=r"$\mathrm{intra} + \mathrm{inter} = 1 - Q_{\max}$")
for name, (color, size) in style.items():
    i = choices[name]
    ax_front.scatter(intra[i], inter[i], s=size, facecolors="none", edgecolors=color,
                     linewidths=2, zorder=3, label=name)
i = choices["oracle (uses labels)"]
ax_front.scatter(intra[i], inter[i], s=110, marker="*", color=INK, zorder=4, label="oracle (uses labels)")
ax_front.set_xlim(-0.03, max(intra) + 0.05)
ax_front.set_ylim(-0.03, 1.05)
ax_front.set_xlabel("intra  (fraction of edges between communities)")
ax_front.set_ylabel("inter  (sum of squared degree shares)")
ax_front.set_title("HP-MOCD's front on the karate club", fontsize=11)
ax_front.legend(loc="upper right", fontsize=8, frameon=False)

names = list(choices)
values = [ami[choices[name]] for name in names]
colors = [style[name][0] for name in names[:-1]] + ["white"]
bars = ax_bar.bar(range(len(names)), values, width=0.6, color=colors, edgecolor=colors[:-1] + [INK],
                  linewidth=1.2, zorder=2)
bars[-1].set_hatch("//")
for x, value in enumerate(values):
    ax_bar.text(x, value + 0.02, f"{value:.3f}", ha="center", va="bottom", fontsize=9, color=INK)
ax_bar.set_xticks(range(len(names)), [name.replace(" (", "\n(").replace(" to ", " to\n") for name in names],
                  fontsize=8)
ax_bar.set_ylim(min(0, min(values)) - 0.02, 1.08)
ax_bar.set_ylabel("AMI against the recorded split")
ax_bar.set_title("Agreement of each choice with the split", fontsize=11)
ax_bar.grid(axis="y", color="#e4e3df", lw=0.8, zorder=0)
for ax in (ax_front, ax_bar):
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
plt.show()

# %%
# The dashed line is a level set of modularity: every point on it has the same
# $Q$ as the maximum-modularity member, and the rest of the front lies on or
# above it. The table lists each choice with its index in the de-duplicated
# front, its number of communities $k$, its modularity and its AMI.

print(f"{'rule':<22} {'member':>6} {'k':>3} {'Q':>7} {'AMI':>6}")
for name, i in choices.items():
    print(f"{name:<22} {i:>6} {k[i]:>3} {Q[i]:>7.4f} {ami[i]:>6.3f}")

# %%
# ## Reading the result
#
# The oracle is the maximum AMI over the front, so no rule applied to this
# front can score above its bar against the same reference; the gap between
# another bar and the oracle's is what that rule leaves behind on this front.
# Several members can share a number of communities with different AMI, so the
# member closest to the target $k$ need not be the recorded split. When two
# rules pick the same member, its rings overlap in the left panel and the two
# rows of the table agree.
#
# One front of one nondeterministic run on a 34-node graph says little about how
# the rules compare in general. [Working with Pareto fronts](../../guide/fronts.md)
# covers each rule in more detail, and
# [RIMPSO's resolution profile](plot_resolution_profile.md) shows the RIMPSO
# selector on RIMPSO's own front.
