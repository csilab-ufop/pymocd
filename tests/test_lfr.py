import math
import unittest
import warnings
from fractions import Fraction

import pymocd

# networkx 3.4.2's messages for the conditions it checks
TAU1 = "tau1 must be greater than one"
TAU2 = "tau2 must be greater than one"
MU = "mu must be in the interval [0, 1]"
MAX_DEGREE = "max_degree must be in the interval (0, n]"
EXACTLY_ONE = "Must assign exactly one of min_degree and average_degree"

N = 2000
BASE = dict(n=N, tau1=2.5, tau2=1.5, min_degree=10, max_degree=40, min_community=40, max_community=120)


def realized(n, edges, gt):
    """Degree, internal degree and community sizes of a generated graph."""
    deg = [0] * n
    kin = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
        if gt[u] == gt[v]:
            kin[u] += 1
            kin[v] += 1
    sizes = {}
    for c in gt.values():
        sizes[c] = sizes.get(c, 0) + 1
    return deg, kin, sizes


def truncated_pmf(tau, lo, hi):
    """p(k) ~ k^-tau on [lo, hi], as a dict."""
    weights = {k: k**-tau for k in range(lo, hi + 1)}
    total = sum(weights.values())
    return {k: w / total for k, w in weights.items()}


def moments(pmf):
    mean = sum(k * p for k, p in pmf.items())
    var = sum((k - mean) ** 2 * p for k, p in pmf.items())
    return mean, math.sqrt(var)


def average_degree_law(tau, average, kmax):
    """The documented law for average_degree: with E(m) the mean of k^-tau on
    [m, kmax] and E(j) <= average < E(j + 1), (1 - a) law[j, kmax] + a law[j + 1, kmax],
    a = (E(j) - average) / (E(j) - E(j + 1))."""
    means = {m: moments(truncated_pmf(tau, m, kmax))[0] for m in range(1, kmax + 1)}
    j = max(m for m, e in means.items() if e <= average)
    a = (means[j] - average) / (means[j] - means[j + 1])
    low, high = truncated_pmf(tau, j, kmax), truncated_pmf(tau, j + 1, kmax)
    return j, {k: (1 - a) * p + a * high.get(k, 0.0) for k, p in low.items()}


class TestLfrGraph(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graphs = {mu: pymocd.lfr_graph(mu=mu, seed=7, **BASE) for mu in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)}

    def test_simple_graph_on_the_requested_nodes(self):
        for mu, (edges, gt) in self.graphs.items():
            with self.subTest(mu=mu):
                self.assertIsInstance(edges, list)
                self.assertIsInstance(gt, dict)
                self.assertEqual(set(gt), set(range(N)))
                self.assertTrue(all(isinstance(u, int) and isinstance(v, int) for u, v in edges))
                self.assertTrue(all(0 <= u < v < N for u, v in edges), "u < v, no self-loop")
                self.assertEqual(edges, sorted(set(edges)), "sorted, every edge once")
                labels = sorted(set(gt.values()))
                self.assertEqual(labels, list(range(len(labels))))

    def test_degrees_within_bounds(self):
        for mu, (edges, gt) in self.graphs.items():
            with self.subTest(mu=mu):
                deg, _, _ = realized(N, edges, gt)
                self.assertGreaterEqual(min(deg), 10)
                self.assertLessEqual(max(deg), 40)

    def test_degree_sequence_does_not_depend_on_mu(self):
        # a draw's degrees and sizes come before anything that depends on mu,
        # and every mu here succeeds on its first draw
        sequences = {mu: realized(N, *g)[0] for mu, g in self.graphs.items()}
        sizes = {mu: sorted(realized(N, *g)[2].values()) for mu, g in self.graphs.items()}
        self.assertEqual(len({tuple(s) for s in sequences.values()}), 1)
        self.assertEqual(len({tuple(s) for s in sizes.values()}), 1)

    def test_community_sizes(self):
        for mu, (edges, gt) in self.graphs.items():
            with self.subTest(mu=mu):
                _, _, sizes = realized(N, edges, gt)
                self.assertEqual(sum(sizes.values()), N)
                self.assertTrue(all(40 <= s <= 120 for s in sizes.values()), sorted(sizes.values()))

    def assert_rounded(self, mu, edges, gt, n):
        """Every internal degree is floor or ceil of the decimal-exact (1 - mu) k;
        an integer target may be missed by one edge, at most once per community.
        Returns the number of such misses."""
        deg, kin, sizes = realized(n, edges, gt)
        keep = 1 - Fraction(str(mu))
        off = 0
        for v in range(n):
            target = keep * deg[v]
            self.assertLess(kin[v], sizes[gt[v]], "internal degree fits the community")
            if target.denominator == 1:
                if kin[v] != target:
                    self.assertEqual(abs(kin[v] - target), 1, (v, deg[v], kin[v]))
                    off += 1
            else:
                self.assertTrue(math.floor(target) <= kin[v] <= math.ceil(target), (v, deg[v], kin[v]))
        self.assertLessEqual(off, len(sizes))
        return off

    def test_internal_degree_is_rounded_mixing(self):
        for mu, (edges, gt) in self.graphs.items():
            with self.subTest(mu=mu):
                self.assertEqual(self.assert_rounded(mu, edges, gt, N), 0)

    def test_integer_targets_are_exact(self):
        # (1 - mu) k is an integer for these pairs, but not in floating point
        # (10 * (1 - 0.7) = 3.0000000000000004): no node may be a flip away
        for mu, k in ((0.7, 10), (0.8, 15), (0.9, 10), (0.3, 90)):
            self.assertNotEqual(k * (1 - mu), round(k * (1 - mu)))
            for seed in (1, 2):
                with self.subTest(mu=mu, seed=seed):
                    edges, gt = pymocd.lfr_graph(1000, 2.5, 1.5, mu, min_degree=k, max_degree=k + 4,
                                                 min_community=k + 10, max_community=k + 80, seed=seed)
                    # every community has fractional members: no stub move
                    self.assertEqual(self.assert_rounded(mu, edges, gt, 1000), 0)

    def test_odd_all_integer_community_moves_one_stub(self):
        # three communities of 33 nodes of degree 10 at mu = 0.5: each needs
        # 33 * 5 internal edge ends, an odd number, so one node per community
        # ends one edge off its target of 5
        edges, gt = pymocd.lfr_graph(99, 2.5, 1.5, 0.5, min_degree=10, max_degree=10, min_community=33,
                                     max_community=33, seed=1)
        deg, kin, sizes = realized(99, edges, gt)
        self.assertEqual(sorted(sizes.values()), [33, 33, 33])
        self.assertEqual(set(deg), {10})
        self.assertEqual(self.assert_rounded(0.5, edges, gt, 99), 3)
        for c in sizes:
            self.assertEqual(sum(kin[v] != 5 for v in range(99) if gt[v] == c), 1)

    def test_realized_mixing_matches_mu(self):
        for mu, (edges, gt) in self.graphs.items():
            with self.subTest(mu=mu):
                deg, kin, _ = realized(N, edges, gt)
                inter = sum(1 for u, v in edges if gt[u] != gt[v])
                self.assertAlmostEqual(inter / len(edges), mu, delta=0.005)
                node_mean = sum((d - i) / d for d, i in zip(deg, kin)) / N
                self.assertAlmostEqual(node_mean, mu, delta=0.005)
                if mu == 0.0:
                    self.assertEqual(inter, 0)
                if mu == 1.0:
                    self.assertEqual(inter, len(edges))

    def test_average_degree_law(self):
        # tau1 = 2.5 on [11, 50] and [12, 50] have means 18.95 and 20.25:
        # the law is the mixture of the two with mean 20
        j, pmf = average_degree_law(2.5, 20, 50)
        self.assertEqual(j, 11)
        mean, sd = moments(pmf)
        self.assertAlmostEqual(mean, 20, delta=1e-9)
        params = dict(n=N, tau1=2.5, tau2=1.5, mu=0.3, average_degree=20, max_degree=50,
                      min_community=20, max_community=100)
        degrees = []
        for seed in range(10):
            edges, gt = pymocd.lfr_graph(**params, seed=seed)
            deg, _, _ = realized(N, edges, gt)
            self.assertEqual(min(deg), j, "p(11) * n = 55 nodes expected at the minimum")
            self.assertLessEqual(max(deg), 50)
            degrees += deg
        # the mean of 20,000 degrees, within 5 standard errors of the request
        self.assertAlmostEqual(sum(degrees) / len(degrees), 20, delta=5 * sd / math.sqrt(len(degrees)))

    def test_low_average_degree_is_met(self):
        # no integer minimum degree comes close: E(1) = 1.80, E(2) = 4.15
        j, pmf = average_degree_law(2.5, 3, 100)
        self.assertEqual(j, 1)
        _, sd = moments(pmf)
        params = dict(n=N, tau1=2.5, tau2=1.5, mu=0.5, average_degree=3, max_degree=100,
                      min_community=60, max_community=200)
        degrees = []
        for seed in range(10):
            edges, gt = pymocd.lfr_graph(**params, seed=seed)
            degrees += realized(N, edges, gt)[0]
        self.assertAlmostEqual(sum(degrees) / len(degrees), 3, delta=5 * sd / math.sqrt(len(degrees)))

    def test_default_community_bounds(self):
        edges, gt = pymocd.lfr_graph(1000, 3.0, 1.5, 0.3, average_degree=20, max_degree=50, seed=2)
        deg, _, sizes = realized(1000, edges, gt)
        self.assertTrue(all(min(deg) <= s <= max(deg) for s in sizes.values()))

    def test_mixing_unbiased_near_zero_and_one(self):
        # communities of 10 to 30 nodes with degrees 3 to 9: rounding the odd
        # member of every odd community the same way would take about half an
        # edge end per community off one side, a z of -20 or so over 40 graphs
        for mu in (0.001, 0.999):
            deviation = variance = 0.0
            keep = 1 - Fraction(str(mu))
            for seed in range(40):
                edges, gt = pymocd.lfr_graph(1000, 2.5, 1.5, mu, min_degree=3, max_degree=9,
                                             min_community=10, max_community=30, seed=seed)
                deg, kin, _ = realized(1000, edges, gt)
                for k, i in zip(deg, kin):
                    target = keep * k
                    f = float(target - math.floor(target))
                    deviation += i - float(target)
                    variance += f * (1 - f)
                variance += 1  # the coin that makes the internal total even
            with self.subTest(mu=mu):
                self.assertLess(abs(deviation / math.sqrt(variance)), 5)

    def test_floor_admission(self):
        # networkx places these nodes on round((1 - mu) k), which is the floor
        # here; the ceiling fits no community
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            edges, gt = pymocd.lfr_graph(2, 2.5, 1.5, 0.5, min_degree=1, seed=1)
        self.assertEqual(edges, [(0, 1)])
        self.assertNotEqual(gt[0], gt[1])
        # degree 9 at mu = 0.5 in communities of 5: 4.5 rounds to 4, as the
        # community only holds 4 internal edges per node. No node can be
        # rounded up for those that rounded up, so the mixing is 5/9, not 0.5,
        # and a warning says so
        with self.assertWarns(RuntimeWarning) as caught:
            edges, gt = pymocd.lfr_graph(60, 2.5, 1.5, 0.5, min_degree=9, max_degree=9,
                                         min_community=5, max_community=5, seed=1)
        deg, kin, sizes = realized(60, edges, gt)
        self.assertEqual(set(sizes.values()), {5})
        self.assertEqual((set(deg), set(kin)), ({9}, {4}))
        self.assertIn(f"graph is {5 / 9:.4f}, against mu=0.5", str(caught.warning))

    def test_mixing_shortfall_warning(self):
        # every degree 2 at mu = 0.3 with the default community bounds, so
        # communities of 2 that hold one internal edge per node: the nodes
        # that rounded 1.4 up are held to 1 and no node can take their edge
        # back, so half the edges leave their community
        cases = [
            (dict(n=10_000, tau1=2.5, tau2=1.5, mu=0.3, min_degree=2, max_degree=2), 0.5),
            # degree 10 in communities of 10 at mu = 0.06: everyone held to 9
            (dict(n=2000, tau1=2.0, tau2=1.5, mu=0.06, min_degree=10, max_degree=10, min_community=10,
                  max_community=10), 0.1),
        ]
        for params, mixing in cases:
            with self.subTest(params=params), warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                edges, gt = pymocd.lfr_graph(**params, seed=3)
                inter = sum(gt[u] != gt[v] for u, v in edges) / len(edges)
                self.assertAlmostEqual(inter, mixing, delta=1e-12)
                self.assertEqual(len(caught), 1, [str(w.message) for w in caught])
                self.assertIs(caught[0].category, RuntimeWarning)
                message = str(caught[0].message)
                self.assertIn(f"graph is {inter:.4f}, against mu={params['mu']}", message)
                self.assertIn(
                    "A larger min_community (with max_community at least as large) avoids this.",
                    message,
                )
        # a request whose communities hold every ceiling does not warn
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            for mu in (0.0, 0.3, 0.6, 1.0):
                pymocd.lfr_graph(mu=mu, seed=5, **BASE)

    def test_same_seed_same_graph(self):
        a = pymocd.lfr_graph(mu=0.4, seed=11, **BASE)
        self.assertEqual(a, pymocd.lfr_graph(mu=0.4, seed=11, **BASE))
        self.assertNotEqual(a[0], pymocd.lfr_graph(mu=0.4, seed=12, **BASE)[0])

    def test_failed_draws_are_drawn_again(self):
        # networkx's docstring example: networkx 3.4.2 returns a graph for 26
        # of seeds 0 to 49 and had not finished after 120 s on 8. About half
        # of the draws cannot be placed (a hub needs a community larger than
        # any drawn); lfr_graph draws again and returns a graph for every seed
        params = dict(n=250, tau1=3, tau2=1.5, mu=0.1, average_degree=5, min_community=20)
        inter = total = 0
        for seed in range(50):
            edges, gt = pymocd.lfr_graph(**params, seed=seed)
            deg, kin, sizes = realized(250, edges, gt)
            self.assertTrue(all(s >= 20 for s in sizes.values()))
            self.assertTrue(all(kin[v] < sizes[gt[v]] for v in range(250)))
            inter += sum(gt[u] != gt[v] for u, v in edges)
            total += len(edges)
            self.assertEqual((edges, gt), pymocd.lfr_graph(**params, seed=seed))
        self.assertAlmostEqual(inter / total, 0.1, delta=0.005)

    def test_dense_external_graphs_are_wired(self):
        # every degree d in two communities of 50: at mu = 1 the external
        # graph is d-regular bipartite, which exists for every d <= 50; at
        # mu = 0.5 the external degree is d / 2. K(33, 33, 33) is the only
        # graph for degree 66 at mu = 1 in three communities of 33
        for n, size, d, mu in ((100, 50, 49, 1.0), (100, 50, 96, 0.5), (99, 33, 66, 1.0)):
            with self.subTest(n=n, size=size, d=d, mu=mu):
                edges, gt = pymocd.lfr_graph(n, 3.0, 1.5, mu, min_degree=d, max_degree=d, min_community=size,
                                             max_community=size, max_iters=1, seed=1)
                deg, kin, sizes = realized(n, edges, gt)
                self.assertEqual(set(sizes.values()), {size})
                self.assertEqual(set(deg), {d})
                self.assertEqual(set(kin), {round(d * (1 - mu))})

    def test_shortfall_warning(self):
        # 41 nodes of degree 5 or 6 (the law for average_degree 5.98 on [5, 6]
        # puts 0.98 on 6) at mu = 0, in communities of 6 and 7: the only sizes
        # summing to 41 are one 6 and five 7s, so only draws with six or more
        # nodes of degree 5 can be placed, about one in 7,000. The graphs
        # returned have a mean degree of 5.85 or less, 2% and 5.8 standard
        # errors (0.022) below the request
        params = dict(n=41, tau1=2.0, tau2=1.5, mu=0.0, average_degree=5.98, max_degree=6, min_community=6,
                      max_community=7, max_iters=100_000)
        for seed in range(3):
            with self.subTest(seed=seed):
                with self.assertWarns(RuntimeWarning) as caught:
                    edges, _ = pymocd.lfr_graph(**params, seed=seed)
                mean = 2 * len(edges) / 41
                message = str(caught.warning)
                self.assertIn(f"{mean:.3f}", message)
                self.assertIn("below average_degree=5.98", message)
                self.assertIn("draws of degrees and community sizes were rejected", message)
                self.assertIn("max_degree", message)
        # requests whose returned graphs are close to the law do not warn
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            for seed in range(10):
                pymocd.lfr_graph(250, 3, 1.5, 0.1, average_degree=5, min_community=20, seed=seed)
                pymocd.lfr_graph(N, 2.5, 1.5, 0.3, average_degree=20, max_degree=50, min_community=20,
                                 max_community=100, seed=seed)

    def test_whole_floats_are_accepted(self):
        # networkx takes floats for these four; an integer value means that integer
        floats = {**BASE, "min_degree": 10.0, "max_degree": 40.0, "min_community": 40.0, "max_community": 120.0}
        self.assertEqual(pymocd.lfr_graph(mu=0.3, seed=3, **floats), pymocd.lfr_graph(mu=0.3, seed=3, **BASE))

    def test_invalid_parameters_raise_value_error(self):
        cases = [
            (dict(tau1=1.0), TAU1),
            (dict(tau1=float("nan")), TAU1),
            (dict(tau2=0.5), TAU2),
            (dict(mu=-0.1), MU),
            (dict(mu=1.1), MU),
            (dict(mu=float("nan")), MU),
            (dict(max_degree=0), MAX_DEGREE),
            (dict(max_degree=-1), MAX_DEGREE),
            (dict(max_degree=N + 1), MAX_DEGREE),
            (dict(average_degree=20), EXACTLY_ONE),
            (dict(min_degree=None), EXACTLY_ONE),
            (dict(min_degree=0), "min_degree must be at least 1"),
            (dict(min_degree=41), "min_degree must not exceed max_degree"),
            (dict(min_community=130), "min_community must not exceed max_community"),
            (dict(max_community=0), "max_community must be at least 1"),
            (dict(n=2**32, max_degree=None), "n must be smaller than 2**32"),
            # values networkx accepts but no power law is defined by
            (dict(tau1=math.inf), "tau1 must be finite"),
            (dict(tau2=math.inf), "tau2 must be finite"),
            (dict(mu=math.inf), MU),
            (dict(max_degree=math.inf), MAX_DEGREE),
            (dict(max_degree=math.nan), MAX_DEGREE),
            (dict(min_degree=math.nan), "min_degree must be finite"),
            (dict(min_degree=math.inf), "min_degree must be finite"),
            (dict(min_degree=None, average_degree=math.nan), "average_degree must be finite"),
            (dict(min_degree=None, average_degree=math.inf), "average_degree must be finite"),
            (dict(tol=math.nan), "tol must be finite"),
            (dict(min_community=math.nan), "min_community must be finite"),
            (dict(max_community=math.inf), "max_community must be finite"),
            (dict(min_degree=10.5), "min_degree must be an integer"),
            (dict(max_degree=39.5), "max_degree must be an integer"),
            (dict(min_community=40.5), "min_community must be an integer"),
        ]
        for change, message in cases:
            params = {**BASE, "mu": 0.3, **change}
            with self.subTest(change=change):
                with self.assertRaises(ValueError) as caught:
                    pymocd.lfr_graph(**params, seed=1)
                self.assertEqual(str(caught.exception), message)

    def test_infeasible_requests_raise_runtime_error(self):
        cases = [
            (dict(min_degree=None, average_degree=45), "Could not match average_degree"),
            (dict(min_community=30, max_community=30, n=100, max_degree=20), "Could not create power law sequence"),
            (dict(min_community=20, max_community=20, n=100, min_degree=25, max_degree=30),
             "Could not assign communities; try increasing min_community"),
            # networkx 3.4.2 raises ExceededMaxIterations with these messages
            (dict(n=100, min_degree=3, max_degree=20, min_community=101, max_community=120),
             "Could not create power law sequence"),
            (dict(n=51, min_degree=51, max_degree=None, min_community=None, max_community=None),
             "Could not create power law sequence"),
            (dict(n=50, min_degree=50, max_degree=None, min_community=None, max_community=None),
             "Could not assign communities; try increasing min_community"),
            # every draw fails: degree 10 at mu = 0 needs a community above 10,
            # and only sizes of 10 sum to 100; raised after max_iters draws
            (dict(n=100, min_degree=10, max_degree=10, min_community=10, max_community=11, max_iters=5),
             "Could not assign communities; try increasing min_community"),
            # a single community cannot host the edges between communities
            (dict(n=60, mu=0.5, min_degree=5, max_degree=10, min_community=60, max_community=60, max_iters=5),
             "Could not wire the inter-community edges; try decreasing max_community"),
        ]
        for change, message in cases:
            params = {**BASE, "mu": 0.0, **change}
            with self.subTest(change=change):
                with self.assertRaises(RuntimeError) as caught:
                    pymocd.lfr_graph(**params, seed=1)
                self.assertEqual(str(caught.exception), message)


    def test_seed_must_fit_an_unsigned_64_bit_integer(self):
        for seed in (-1, 2**64):
            with self.subTest(seed=seed), self.assertRaises(OverflowError):
                pymocd.lfr_graph(mu=0.3, seed=seed, **BASE)


if __name__ == "__main__":
    unittest.main()
