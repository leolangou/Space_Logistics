import sys
import unittest
from pathlib import Path

PACKET_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKET_ROOT))

from pipeline import (  # noqa: E402
    ObjectiveWeights,
    OrbitSite,
    graph_betti_numbers,
    optimize_placement,
    transfer_cost_km_s,
)


class TransferCostTests(unittest.TestCase):
    def test_identical_circular_orbits_have_zero_surrogate_cost(self):
        orbit = OrbitSite("same", 550.0, 53.0, 20.0)
        self.assertAlmostEqual(
            transfer_cost_km_s(orbit, orbit, round_trip=False),
            0.0,
            places=12,
        )

    def test_plane_change_cost_is_symmetric_and_positive(self):
        first = OrbitSite("first", 550.0, 0.0, 0.0)
        second = OrbitSite("second", 550.0, 30.0, 0.0)
        forward = transfer_cost_km_s(first, second, round_trip=False)
        reverse = transfer_cost_km_s(second, first, round_trip=False)
        self.assertGreater(forward, 0.0)
        self.assertAlmostEqual(forward, reverse, places=12)


class TopologyTests(unittest.TestCase):
    def test_bipartite_cycle_appears_with_redundant_service(self):
        costs = [[1.0, 3.0], [2.0, 1.0]]
        sparse = graph_betti_numbers(costs, (0, 1), 1.5)
        dense = graph_betti_numbers(costs, (0, 1), 3.0)
        self.assertEqual((sparse["beta_0"], sparse["beta_1"]), (2, 0))
        self.assertEqual((dense["beta_0"], dense["beta_1"]), (1, 1))

    def test_optimizer_returns_requested_number_of_depots(self):
        customers = [
            OrbitSite("C1", 0.0, 0.0, 0.0),
            OrbitSite("C2", 0.0, 0.0, 0.0),
        ]
        candidates = [
            OrbitSite("D1", 0.0, 0.0, 0.0),
            OrbitSite("D2", 0.0, 0.0, 0.0),
            OrbitSite("D3", 0.0, 0.0, 0.0),
        ]
        costs = [[1.0, 2.0, 9.0], [2.0, 1.0, 9.0]]
        selected, metrics = optimize_placement(
            costs,
            customers,
            candidates,
            depot_count=2,
            redundancy=2,
            coverage_threshold_km_s=3.0,
            weights=ObjectiveWeights(),
        )
        self.assertEqual(selected, (0, 1))
        self.assertEqual(len(selected), 2)
        self.assertEqual(metrics["weighted_coverage_deficit"], 0.0)


if __name__ == "__main__":
    unittest.main()
