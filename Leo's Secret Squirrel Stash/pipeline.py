"""Cadet-scale baseline for TDA-informed orbital refueling-depot placement.

The transfer model is intentionally low fidelity. Replace transfer_cost_km_s()
with a validated Lambert, low-thrust, or trajectory-optimization model before
drawing operational conclusions.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

MU_EARTH_KM3_S2 = 398600.4418
R_EARTH_KM = 6378.1363


@dataclass(frozen=True)
class OrbitSite:
    name: str
    altitude_km: float
    inclination_deg: float
    raan_deg: float
    weight: float = 1.0


@dataclass(frozen=True)
class ObjectiveWeights:
    backup: float = 0.35
    single_failure: float = 0.25
    coverage_deficit: float = 20.0
    disconnected_component: float = 1.0
    cycle_reward: float = 0.05


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size
        self.components = size

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> bool:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return False
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1
        self.components -= 1
        return True


def load_sites(path: Path) -> list[OrbitSite]:
    """Load orbit sites from a CSV with name, altitude, inclination, and RAAN."""
    sites: list[OrbitSite] = []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            sites.append(
                OrbitSite(
                    name=row["name"],
                    altitude_km=float(row["altitude_km"]),
                    inclination_deg=float(row["inclination_deg"]),
                    raan_deg=float(row["raan_deg"]),
                    weight=float(row.get("weight") or 1.0),
                )
            )
    if not sites:
        raise ValueError(f"No orbit sites found in {path}")
    return sites


def relative_plane_angle_rad(first: OrbitSite, second: OrbitSite) -> float:
    """Angle between two orbit-plane normals."""
    i1 = math.radians(first.inclination_deg)
    i2 = math.radians(second.inclination_deg)
    delta_raan = math.radians(first.raan_deg - second.raan_deg)
    cosine = (
        math.cos(i1) * math.cos(i2)
        + math.sin(i1) * math.sin(i2) * math.cos(delta_raan)
    )
    return math.acos(max(-1.0, min(1.0, cosine)))


def transfer_cost_km_s(
    first: OrbitSite,
    second: OrbitSite,
    *,
    round_trip: bool = True,
) -> float:
    """Return a transparent circular-orbit transfer-cost surrogate.

    The estimate adds a coplanar Hohmann transfer to a plane change performed at
    the larger orbital radius. It omits phasing, RAAN drift, rendezvous timing,
    finite-burn effects, perturbations, and propellant logistics.
    """
    radius_1 = R_EARTH_KM + first.altitude_km
    radius_2 = R_EARTH_KM + second.altitude_km
    transfer_semimajor_axis = 0.5 * (radius_1 + radius_2)

    circular_1 = math.sqrt(MU_EARTH_KM3_S2 / radius_1)
    circular_2 = math.sqrt(MU_EARTH_KM3_S2 / radius_2)
    transfer_1 = math.sqrt(
        MU_EARTH_KM3_S2 * (2.0 / radius_1 - 1.0 / transfer_semimajor_axis)
    )
    transfer_2 = math.sqrt(
        MU_EARTH_KM3_S2 * (2.0 / radius_2 - 1.0 / transfer_semimajor_axis)
    )
    hohmann = abs(transfer_1 - circular_1) + abs(circular_2 - transfer_2)

    plane_angle = relative_plane_angle_rad(first, second)
    plane_change_speed = math.sqrt(
        MU_EARTH_KM3_S2 / max(radius_1, radius_2)
    )
    plane_change = 2.0 * plane_change_speed * math.sin(0.5 * plane_angle)

    one_way = hohmann + plane_change
    return 2.0 * one_way if round_trip else one_way


def build_cost_matrix(
    customers: Sequence[OrbitSite],
    candidates: Sequence[OrbitSite],
    *,
    round_trip: bool = True,
) -> list[list[float]]:
    return [
        [
            transfer_cost_km_s(customer, candidate, round_trip=round_trip)
            for candidate in candidates
        ]
        for customer in customers
    ]


def graph_betti_numbers(
    cost_matrix: Sequence[Sequence[float]],
    selected: Sequence[int],
    epsilon_km_s: float,
) -> dict[str, int]:
    """Betti numbers of the thresholded bipartite service graph.

    Customer and selected-depot vertices are joined when transfer cost is at
    most epsilon. For a graph, beta_0 is the component count and
    beta_1 = edges - vertices + components is the independent-cycle count.
    """
    customer_count = len(cost_matrix)
    depot_count = len(selected)
    vertex_count = customer_count + depot_count
    union_find = UnionFind(vertex_count)
    edge_count = 0

    for customer_index, row in enumerate(cost_matrix):
        for local_depot_index, candidate_index in enumerate(selected):
            if row[candidate_index] <= epsilon_km_s:
                edge_count += 1
                union_find.union(
                    customer_index,
                    customer_count + local_depot_index,
                )

    beta_0 = union_find.components
    beta_1 = edge_count - vertex_count + beta_0
    return {
        "epsilon_km_s": epsilon_km_s,
        "vertices": vertex_count,
        "edges": edge_count,
        "beta_0": beta_0,
        "beta_1": beta_1,
    }


def betti_curve(
    cost_matrix: Sequence[Sequence[float]],
    selected: Sequence[int],
    epsilons_km_s: Iterable[float],
) -> list[dict[str, int]]:
    return [
        graph_betti_numbers(cost_matrix, selected, epsilon)
        for epsilon in epsilons_km_s
    ]


def weighted_mean(values: Sequence[float], customers: Sequence[OrbitSite]) -> float:
    total_weight = sum(customer.weight for customer in customers)
    return sum(
        customer.weight * value
        for customer, value in zip(customers, values)
    ) / total_weight


def evaluate_placement(
    cost_matrix: Sequence[Sequence[float]],
    customers: Sequence[OrbitSite],
    selected: Sequence[int],
    *,
    redundancy: int,
    coverage_threshold_km_s: float,
    weights: ObjectiveWeights,
) -> dict[str, object]:
    if redundancy < 1:
        raise ValueError("redundancy must be at least 1")
    if len(selected) < redundancy:
        raise ValueError("selected depot count must be at least redundancy")

    selected_costs = [
        sorted(row[index] for index in selected)
        for row in cost_matrix
    ]
    primary_costs = [row[0] for row in selected_costs]
    backup_costs = [row[redundancy - 1] for row in selected_costs]
    coverage_counts = [
        sum(cost <= coverage_threshold_km_s for cost in row)
        for row in selected_costs
    ]
    deficits = [
        max(0, redundancy - count)
        for count in coverage_counts
    ]

    failure_means: list[float] = []
    for failed in selected:
        surviving_primary = [
            min(row[index] for index in selected if index != failed)
            for row in cost_matrix
        ]
        failure_means.append(weighted_mean(surviving_primary, customers))

    topology = graph_betti_numbers(
        cost_matrix,
        selected,
        coverage_threshold_km_s,
    )
    primary_mean = weighted_mean(primary_costs, customers)
    backup_mean = weighted_mean(backup_costs, customers)
    failure_worst = max(failure_means)
    deficit_mean = weighted_mean(deficits, customers)

    objective = (
        primary_mean
        + weights.backup * backup_mean
        + weights.single_failure * failure_worst
        + weights.coverage_deficit * deficit_mean
        + weights.disconnected_component * max(0, topology["beta_0"] - 1)
        - weights.cycle_reward * topology["beta_1"]
    )

    return {
        "objective": objective,
        "mean_primary_delta_v_km_s": primary_mean,
        "mean_redundant_delta_v_km_s": backup_mean,
        "worst_single_failure_mean_delta_v_km_s": failure_worst,
        "weighted_coverage_deficit": deficit_mean,
        "coverage_counts": coverage_counts,
        "topology_at_threshold": topology,
    }


def optimize_placement(
    cost_matrix: Sequence[Sequence[float]],
    customers: Sequence[OrbitSite],
    candidates: Sequence[OrbitSite],
    *,
    depot_count: int,
    redundancy: int,
    coverage_threshold_km_s: float,
    weights: ObjectiveWeights,
) -> tuple[tuple[int, ...], dict[str, object]]:
    if depot_count > len(candidates):
        raise ValueError("depot_count exceeds the number of candidates")
    if depot_count < redundancy:
        raise ValueError("depot_count must be at least redundancy")

    best_selection: tuple[int, ...] | None = None
    best_metrics: dict[str, object] | None = None
    for selected in itertools.combinations(range(len(candidates)), depot_count):
        metrics = evaluate_placement(
            cost_matrix,
            customers,
            selected,
            redundancy=redundancy,
            coverage_threshold_km_s=coverage_threshold_km_s,
            weights=weights,
        )
        if best_metrics is None or metrics["objective"] < best_metrics["objective"]:
            best_selection = selected
            best_metrics = metrics

    assert best_selection is not None and best_metrics is not None
    return best_selection, best_metrics


def evenly_spaced(start: float, stop: float, count: int) -> list[float]:
    if count < 2:
        return [start]
    step = (stop - start) / (count - 1)
    return [start + index * step for index in range(count)]


def main() -> None:
    packet_root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Optimize static orbital depot placement and compute Betti curves."
    )
    parser.add_argument(
        "--customers",
        type=Path,
        default=packet_root / "data" / "example_customers.csv",
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=packet_root / "data" / "example_candidates.csv",
    )
    parser.add_argument("--depots", type=int, default=3)
    parser.add_argument("--redundancy", type=int, default=2)
    parser.add_argument("--coverage-threshold", type=float, default=12.0)
    parser.add_argument("--one-way", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    customers = load_sites(args.customers)
    candidates = load_sites(args.candidates)
    matrix = build_cost_matrix(
        customers,
        candidates,
        round_trip=not args.one_way,
    )
    weights = ObjectiveWeights()
    selected, metrics = optimize_placement(
        matrix,
        customers,
        candidates,
        depot_count=args.depots,
        redundancy=args.redundancy,
        coverage_threshold_km_s=args.coverage_threshold,
        weights=weights,
    )

    curve = betti_curve(
        matrix,
        selected,
        evenly_spaced(0.0, 20.0, 41),
    )
    result = {
        "model_warning": (
            "Educational surrogate only; replace the transfer-cost model before "
            "using results for mission design."
        ),
        "selected_depots": [candidates[index].name for index in selected],
        "configuration": {
            "depot_count": args.depots,
            "redundancy": args.redundancy,
            "coverage_threshold_km_s": args.coverage_threshold,
            "round_trip": not args.one_way,
            "objective_weights": asdict(weights),
        },
        "metrics": metrics,
        "betti_curve": curve,
    }
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
