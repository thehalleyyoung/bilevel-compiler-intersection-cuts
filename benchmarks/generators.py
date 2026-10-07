"""Random bilevel instance generators used by kkt_eval.py.

Each generator returns a BilevelInstance

    min_{x,y}  d^T x + e^T y
    s.t.       C x + D y <= h,  x_lb <= x <= x_ub,  x integer
               y in argmin { c^T y' : A y' <= b + B x, y' >= 0 }

with a continuous (LP) follower. The four families are a knapsack-type
problem with an interdicted capacity, dense random bilevel programs,
programs with pairwise follower coupling, and a Stackelberg
target-defence game with an LP attacker.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class BilevelInstance:
    """Bilevel instance with integer leader and LP follower."""
    name: str
    n_x: int
    n_y: int
    m_upper: int
    m_lower: int
    d: np.ndarray
    e: np.ndarray
    C: np.ndarray
    D: np.ndarray
    h: np.ndarray
    c: np.ndarray
    A: np.ndarray
    b: np.ndarray
    B: np.ndarray
    x_lb: np.ndarray = None
    x_ub: np.ndarray = None
    x_binary: bool = True
    category: str = "general"

    def __post_init__(self):
        if self.x_lb is None:
            self.x_lb = np.zeros(self.n_x)
        if self.x_ub is None:
            self.x_ub = np.ones(self.n_x)


def generate_hard_knapsack_interdiction(n: int, seed: int = 42,
                                         correlation: str = "weakly") -> BilevelInstance:
    """Knapsack-type instance: binary x removes capacity from an LP knapsack.

    Weights and values are weakly correlated, strongly correlated or
    uncorrelated (the three classes of Pisinger's knapsack generator).
    Leader and follower share the objective -values^T y; the leader has a
    budget constraint and a joint constraint sum(x) + sum(y) <= 0.7 n.
    """
    rng = np.random.RandomState(seed)

    if correlation == "weakly":
        weights = rng.randint(10, 100, size=n).astype(float)
        values = weights + rng.randint(-10, 11, size=n).astype(float)
        values = np.maximum(values, 1.0)
    elif correlation == "strongly":
        weights = rng.randint(10, 100, size=n).astype(float)
        values = weights + 10.0
    else:
        weights = rng.randint(1, 50, size=n).astype(float)
        values = rng.randint(1, 50, size=n).astype(float)

    # Tight capacity → harder
    capacity = 0.35 * weights.sum()
    budget = max(1, n // 4)

    n_x = n
    n_y = n
    m_upper = 2  # budget + linking
    m_lower = 1 + n  # capacity + upper bounds

    # Leader: min -values^T y (maximize interdiction impact)
    d = np.zeros(n_x)
    e = -values.copy()

    # Upper-level: budget + a linking constraint to create coupling
    C_mat = np.zeros((m_upper, n_x))
    D_mat = np.zeros((m_upper, n_y))
    h_vec = np.zeros(m_upper)
    C_mat[0, :] = 1.0
    h_vec[0] = budget
    # Linking: sum(x_i + y_i) <= n * 0.7 (creates interdependence)
    C_mat[1, :] = 1.0
    D_mat[1, :] = 1.0
    h_vec[1] = n * 0.7

    c_vec = -values.copy()

    A_mat = np.zeros((m_lower, n_y))
    b_vec = np.zeros(m_lower)
    B_mat = np.zeros((m_lower, n_x))

    # Capacity with interdiction
    A_mat[0, :] = weights
    b_vec[0] = capacity
    B_mat[0, :] = -weights  # capacity shrinks with interdiction

    # y_i <= 1
    for i in range(n):
        A_mat[1 + i, i] = 1.0
        b_vec[1 + i] = 1.0

    return BilevelInstance(
        name=f"hard_knap_n{n}_{correlation[:4]}_s{seed}",
        n_x=n_x, n_y=n_y, m_upper=m_upper, m_lower=m_lower,
        d=d, e=e, C=C_mat, D=D_mat, h=h_vec,
        c=c_vec, A=A_mat, b=b_vec, B=B_mat,
        x_binary=True, category="knapsack_interdiction",
    )


def generate_dense_bilevel(n_x: int, n_y: int, density: float = 0.7,
                            seed: int = 42) -> BilevelInstance:
    """Dense random bilevel program; x is integer in [0, 3]."""
    rng = np.random.RandomState(seed)

    m_upper = max(2, n_x // 2)
    m_lower = max(3, n_y + n_x // 2)

    d = rng.uniform(-5, 5, size=n_x)
    e = rng.uniform(-5, 5, size=n_y)

    # Dense upper-level constraints
    C = rng.uniform(-1, 3, size=(m_upper, n_x))
    C[rng.random(C.shape) > density] = 0
    D = rng.uniform(-1, 3, size=(m_upper, n_y))
    D[rng.random(D.shape) > density] = 0
    h = rng.uniform(5, 20, size=m_upper)

    c = rng.uniform(-3, 3, size=n_y)

    # Dense lower-level
    A = rng.uniform(0, 3, size=(m_lower, n_y))
    A[rng.random(A.shape) > density] = 0
    b = rng.uniform(2, 15, size=m_lower)
    B = rng.uniform(-2, 2, size=(m_lower, n_x))
    B[rng.random(B.shape) > density] = 0

    return BilevelInstance(
        name=f"dense_{n_x}x{n_y}_d{int(density*10)}_s{seed}",
        n_x=n_x, n_y=n_y, m_upper=m_upper, m_lower=m_lower,
        d=d, e=e, C=C, D=D, h=h,
        c=c, A=A, b=b, B=B,
        x_binary=True,
        x_ub=np.full(n_x, 3.0),
        category="dense_bilevel",
    )


def generate_bilevel_with_integer_linking(n: int, seed: int = 42) -> BilevelInstance:
    """Binary leader, LP follower with pairwise coupling y_i + 0.5 y_{i+1}."""
    rng = np.random.RandomState(seed)

    n_x = n
    n_y = n
    m_upper = 2
    m_lower = 2 * n  # many constraints to create more duals

    d = rng.uniform(-3, 1, size=n_x)
    e = rng.uniform(-5, -1, size=n_y)  # leader benefits from follower activity

    C_mat = np.zeros((m_upper, n_x))
    D_mat = np.zeros((m_upper, n_y))
    h_vec = np.zeros(m_upper)
    C_mat[0, :] = 1.0
    h_vec[0] = n * 0.6
    D_mat[1, :] = 1.0
    h_vec[1] = n * 0.8

    c_vec = rng.uniform(1, 10, size=n_y)

    A_mat = np.zeros((m_lower, n_y))
    b_vec = np.zeros(m_lower)
    B_mat = np.zeros((m_lower, n_x))

    # Pairwise coupling constraints: y_i + y_{i+1} <= b_k + B_k x
    for k in range(n):
        A_mat[k, k] = 1.0
        if k + 1 < n:
            A_mat[k, k + 1] = 0.5
        b_vec[k] = rng.uniform(1, 3)
        B_mat[k, k] = -rng.uniform(0.2, 0.8)

    # Individual bounds
    for k in range(n):
        A_mat[n + k, k] = 1.0
        b_vec[n + k] = rng.uniform(2, 5)
        B_mat[n + k, min(k, n_x - 1)] = -rng.uniform(0.1, 0.5)

    return BilevelInstance(
        name=f"intlink_n{n}_s{seed}",
        n_x=n_x, n_y=n_y, m_upper=m_upper, m_lower=m_lower,
        d=d, e=e, C=C_mat, D=D_mat, h=h_vec,
        c=c_vec, A=A_mat, b=b_vec, B=B_mat,
        x_binary=True, category="integer_linking",
    )


def generate_stackelberg_game(n_strategies: int, seed: int = 42) -> BilevelInstance:
    """Target-defence game: binary x defends targets, an LP attacker picks y.

    Defending target i lowers its attack bound from 1 to 0.5; the attacker
    maximises damage^T y subject to sum(y) <= 1 and the leader minimises it.
    """
    rng = np.random.RandomState(seed)

    n_targets = n_strategies
    n_x = n_targets  # which targets to defend
    n_y = n_targets  # which target to attack (probability)

    # Leader defends at most k targets
    k = max(1, n_targets // 3)

    m_upper = 1
    m_lower = n_targets + 1  # attack probabilities + normalization

    # Leader: min sum(attack_damage * (1 - x_i) * y_i)
    # Linearized: min sum(-damage_i * y_i) + sum(damage_i * x_i * y_i)
    # With continuous follower: we use bilinear McCormick or simplified version
    damage = rng.uniform(5, 50, size=n_targets)
    d = np.zeros(n_x)  # leader cost
    e = damage.copy()  # attack damage

    C_mat = np.ones((1, n_x))
    D_mat = np.zeros((1, n_y))
    h_vec = np.array([k])

    # Follower: max damage * (1-x) * y = min -damage * (1-x) * y
    # Simplified (for LP follower): min -damage^T y (attacker maximizes expected damage)
    c_vec = -damage.copy()

    # Follower constraints: probability simplex + capacity
    A_mat = np.zeros((m_lower, n_y))
    b_vec = np.zeros(m_lower)
    B_mat = np.zeros((m_lower, n_x))

    # y_i <= 1 for each target
    for i in range(n_targets):
        A_mat[i, i] = 1.0
        b_vec[i] = 1.0
        B_mat[i, i] = -0.5  # defending reduces attack effectiveness

    # Sum of y <= 1 (probability constraint)
    A_mat[n_targets, :] = 1.0
    b_vec[n_targets] = 1.0

    return BilevelInstance(
        name=f"stackelberg_n{n_targets}_s{seed}",
        n_x=n_x, n_y=n_y, m_upper=m_upper, m_lower=m_lower,
        d=d, e=e, C=C_mat, D=D_mat, h=h_vec,
        c=c_vec, A=A_mat, b=b_vec, B=B_mat,
        x_binary=True, category="stackelberg_game",
    )
