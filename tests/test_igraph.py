import unittest

import networkx as nx
import pymocd

try:
    import igraph as ig
except ImportError:
    ig = None

DETECTORS = [
    "rimpso",
    "hpmocd",
    "cdrme",
    "mmcomo",
    "ccm",
    "krm",
    "gdpso",
    "mocd_q",
    "mocd_d",
    "moga_net",
]


@unittest.skipIf(ig is None, "python-igraph is not installed")
class TestIgraphInput(unittest.TestCase):
    def setUp(self):
        self.graph = ig.Graph.Famous("Zachary")

    def test_every_detector_accepts_igraph(self):
        for name in DETECTORS:
            with self.subTest(detector=name):
                result = getattr(pymocd, name)(self.graph)
                self.assertEqual(set(result), set(range(self.graph.vcount())))

    def test_rimpso_matches_networkx_with_same_order(self):
        twin = nx.Graph()
        twin.add_nodes_from(range(self.graph.vcount()))
        twin.add_edges_from(self.graph.get_edgelist())
        self.assertEqual(pymocd.rimpso(self.graph), pymocd.rimpso(twin))

    def test_isolated_vertex_gets_minus_one(self):
        self.graph.add_vertices(1)
        self.assertEqual(pymocd.rimpso(self.graph)[34], -1)


if __name__ == "__main__":
    unittest.main()
