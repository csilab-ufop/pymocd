import random
import unittest
import warnings

import pymocd

try:
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        normalized_mutual_info_score,
    )
except ImportError:
    adjusted_mutual_info_score = None

TOL = 1e-10


def as_dict(labels):
    return dict(enumerate(labels))


def random_labels(rng, n, k):
    return [rng.randrange(k) for _ in range(n)]


@unittest.skipIf(adjusted_mutual_info_score is None, "scikit-learn is not installed")
class TestMetricsMatchSklearn(unittest.TestCase):
    def assert_matches(self, partition, gt):
        # sklearn scores the shared nodes; labels_true = gt, labels_pred = partition
        nodes = [v for v in partition if v in gt]
        pred = [partition[v] for v in nodes]
        true = [gt[v] for v in nodes]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # many-class labelings warn about regression targets
            expected = (
                normalized_mutual_info_score(true, pred),
                adjusted_mutual_info_score(true, pred),
                adjusted_rand_score(true, pred),
            )
        got = (pymocd.nmi(partition, gt), pymocd.ami(partition, gt), pymocd.ari(partition, gt))
        for name, mine, theirs in zip(("nmi", "ami", "ari"), got, expected):
            self.assertAlmostEqual(mine, theirs, delta=TOL, msg=name)
        self.assertEqual(pymocd.ari(partition, gt), expected[2])
        self.assertEqual(pymocd.gt_metrics(partition, gt), got + (pymocd.f1(partition, gt),))

    def test_random_labelings(self):
        # stops at n = 1000: at n = 2000 with 2000 labels drawn per side,
        # sklearn's own AMI is 2e-10 from the exact value (pymocd's is 2e-16)
        rng = random.Random(7)
        for n in (2, 3, 10, 57, 300, 1000):
            for k_true in (1, 2, 5, n // 3 + 1, n):
                for k_pred in (1, 3, n // 2 + 1, n):
                    with self.subTest(n=n, k_true=k_true, k_pred=k_pred):
                        gt = as_dict(random_labels(rng, n, k_true))
                        partition = as_dict(random_labels(rng, n, k_pred))
                        self.assert_matches(partition, gt)

    def test_noisy_copies(self):
        rng = random.Random(11)
        for n, k in ((200, 4), (1000, 30), (3000, 300)):
            truth = random_labels(rng, n, k)
            for noise in (0.0, 0.1, 0.5):
                with self.subTest(n=n, k=k, noise=noise):
                    pred = [rng.randrange(k) if rng.random() < noise else c for c in truth]
                    self.assert_matches(as_dict(pred), as_dict(truth))

    def test_identical_singletons(self):
        # sklearn returns 1.0; earlier versions of pymocd returned 0.0 at n = 2 and 4
        for n in (2, 4, 14, 22, 309):
            with self.subTest(n=n):
                singletons = as_dict(range(n))
                self.assertEqual(pymocd.ami(singletons, dict(singletons)), 1.0)
                self.assert_matches(singletons, as_dict(reversed(range(n))))

    def test_degenerate_labelings(self):
        cases = {
            "one node": ({0: 0}, {0: 7}),
            "two nodes": ({0: 0, 1: 1}, {0: 0, 1: 1}),
            "both one community": ({v: 0 for v in range(5)}, {v: 4 for v in range(5)}),
            "one community vs singletons": ({v: 0 for v in range(6)}, as_dict(range(6))),
            "one side all singletons": (as_dict(range(6)), as_dict([0, 1, 2, 3, 2, 4])),
            "giant community": (as_dict([0] * 300 + [1]), as_dict([0] * 299 + [1, 1])),
        }
        for name, (partition, gt) in cases.items():
            with self.subTest(case=name):
                self.assert_matches(partition, gt)
                self.assert_matches(gt, partition)

    def test_only_shared_nodes_are_scored(self):
        rng = random.Random(3)
        partition = {v: rng.randrange(4) for v in range(0, 300)}
        gt = {v: rng.randrange(6) for v in range(100, 500)}
        self.assert_matches(partition, gt)

    def test_any_hashable_label_or_node(self):
        rng = random.Random(5)
        big = {v: 2**40 + rng.randrange(4) for v in range(200)}  # beyond 32 bits
        names = {v: rng.choice("abcde") for v in range(200)}
        self.assert_matches(big, {v: rng.randrange(3) for v in range(200)})
        self.assert_matches(names, big)
        self.assert_matches({f"n{v}": c for v, c in names.items()},
                            {f"n{v}": float(c % 3) for v, c in big.items()})


class Index:
    """Converts to an int through __index__ but hashes and compares by identity."""

    def __init__(self, value):
        self.value = value

    def __index__(self):
        return self.value


FUNCTIONS = (pymocd.gt_metrics, pymocd.nmi, pymocd.ami, pymocd.ari, pymocd.f1)


class TestMetricConventions(unittest.TestCase):
    def test_two_empty_dicts_score_one(self):
        single = tuple(f({}, {}) for f in FUNCTIONS[1:])
        self.assertEqual(single, (1.0, 1.0, 1.0, 1.0))
        self.assertEqual(pymocd.gt_metrics({}, {}), single)

    def test_no_shared_node_raises(self):
        cases = {
            "disjoint ids": ({0: 0, 1: 1}, {5: 0, 6: 1}),
            "int ids against str ids": ({0: 0, 1: 0, 2: 1}, {"0": 0, "1": 0, "2": 1}),
            "empty partition": ({}, {0: 1}),
            "empty gt": ({0: 1}, {}),
        }
        for name, (partition, gt) in cases.items():
            for f in FUNCTIONS:
                with self.subTest(case=name, function=f.__name__):
                    with self.assertRaises(ValueError):
                        f(partition, gt)

    def test_one_shared_node_scores_one(self):
        self.assertEqual(pymocd.gt_metrics({0: 0, 1: 1}, {1: 5, 2: 6}), (1.0, 1.0, 1.0, 1.0))

    def test_index_objects_are_matched_as_dict_keys(self):
        # Index(1) != 1 as a dict key, so each Index label is its own community
        # and an Index node id matches no int node id
        gt = {v: v % 2 for v in range(8)}
        partition = {v: Index(v % 2) for v in range(8)}
        as_singletons = {v: v for v in range(8)}
        self.assertEqual(pymocd.gt_metrics(partition, gt), pymocd.gt_metrics(as_singletons, gt))
        # an entry that is not scored changes nothing
        self.assertEqual(pymocd.gt_metrics(partition | {100: "unrelated"}, gt),
                         pymocd.gt_metrics(partition, gt))
        with self.assertRaises(ValueError):
            pymocd.gt_metrics({Index(v): v % 3 for v in range(6)}, {v: v % 2 for v in range(6)})

    def test_numpy_integers_match_ints(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not installed")
        rng = random.Random(9)
        gt = {v: rng.randrange(5) for v in range(300)}
        partition = {v: rng.randrange(7) for v in range(300)}
        as_numpy = {np.int64(v): np.uint8(c) for v, c in partition.items()}
        self.assertEqual(pymocd.gt_metrics(as_numpy, gt), pymocd.gt_metrics(partition, gt))

    def test_minus_one_is_an_ordinary_label(self):
        partition = {0: -1, 1: -1, 2: 0, 3: 0}
        gt = {0: 5, 1: 5, 2: 6, 3: 6}
        self.assertEqual(pymocd.gt_metrics(partition, gt), (1.0, 1.0, 1.0, 1.0))

    def test_equal_labels_name_one_community(self):
        # 1, 1.0 and True compare equal, as dict keys do
        partition = {0: 1, 1: 1.0, 2: True, 3: 0}
        gt = {0: "a", 1: "a", 2: "a", 3: "b"}
        self.assertEqual(pymocd.gt_metrics(partition, gt), (1.0, 1.0, 1.0, 1.0))

    def test_f1_is_one_without_co_clustered_pairs(self):
        self.assertEqual(pymocd.f1({0: 0, 1: 1, 2: 2}, {0: 5, 1: 6, 2: 7}), 1.0)
        self.assertEqual(pymocd.f1({0: 0}, {0: 1}), 1.0)

    def test_unhashable_label_raises(self):
        with self.assertRaises(TypeError):
            pymocd.nmi({0: [1]}, {0: 0})


if __name__ == "__main__":
    unittest.main()
