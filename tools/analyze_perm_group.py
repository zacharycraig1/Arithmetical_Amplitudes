import re
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Cycle = Tuple[int, ...]


def parse_gap_cycles(s: str, n: int) -> List[int]:
    """
    Parse a GAP-style permutation like "(1,2)(3,4,5)" into 1-indexed array perm[1..n].
    Returns a list of length n+1 where index 0 is unused.
    """
    s = s.strip()
    if s in ("()", ""):
        return list(range(n + 1))
    # remove trailing semicolon if present
    if s.endswith(";"):
        s = s[:-1]
    # Extract cycles
    cycles = re.findall(r"\(([^()]*)\)", s)
    perm = list(range(n + 1))
    for cyc in cycles:
        cyc = cyc.strip()
        if not cyc:
            continue
        nums = [int(x.strip()) for x in cyc.split(",") if x.strip()]
        if len(nums) < 2:
            continue
        for a, b in zip(nums, nums[1:]):
            perm[a] = b
        perm[nums[-1]] = nums[0]
    return perm


def compose(p: Sequence[int], q: Sequence[int]) -> List[int]:
    """Return composition p ∘ q (apply q then p), 1-indexed."""
    n = len(p) - 1
    return [0] + [p[q[i]] for i in range(1, n + 1)]


def inv(p: Sequence[int]) -> List[int]:
    n = len(p) - 1
    out = [0] + [0] * n
    for i in range(1, n + 1):
        out[p[i]] = i
    return out


def cycle_decomposition(p: Sequence[int]) -> List[Cycle]:
    n = len(p) - 1
    seen = [False] * (n + 1)
    out: List[Cycle] = []
    for i in range(1, n + 1):
        if seen[i]:
            continue
        if p[i] == i:
            seen[i] = True
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = p[j]
        if len(cyc) > 1:
            out.append(tuple(cyc))
    return out


def perm_order(p: Sequence[int]) -> int:
    from math import gcd

    def lcm(a: int, b: int) -> int:
        return a // gcd(a, b) * b

    ord_ = 1
    for cyc in cycle_decomposition(p):
        ord_ = lcm(ord_, len(cyc))
    return ord_


def cycle_type(p: Sequence[int]) -> Tuple[int, ...]:
    n = len(p) - 1
    cycs = cycle_decomposition(p)
    moved = sum(len(c) for c in cycs)
    fixed = n - moved
    lens = sorted([len(c) for c in cycs], reverse=True)
    if fixed:
        lens.extend([1] * fixed)
    return tuple(sorted(lens, reverse=True))


def block_decomposition_from_tau(tau: Sequence[int]) -> List[Tuple[int, int]]:
    """Return 12 blocks (pairs) from a fixed-point-free involution tau."""
    n = len(tau) - 1
    seen = set()
    blocks = []
    for i in range(1, n + 1):
        if i in seen:
            continue
        j = tau[i]
        if j == i:
            raise ValueError("tau has a fixed point; cannot form 2-blocks")
        seen.add(i)
        seen.add(j)
        blocks.append((i, j))
    return blocks


def induced_block_permutation(perm: Sequence[int], blocks: List[Tuple[int, int]]) -> List[int]:
    """
    Induced permutation on blocks. Blocks are indexed 1..m in order given.
    """
    block_index: Dict[int, int] = {}
    for idx, (a, b) in enumerate(blocks, start=1):
        block_index[a] = idx
        block_index[b] = idx
    m = len(blocks)
    out = [0] + [0] * m
    for idx, (a, _) in enumerate(blocks, start=1):
        image_block = block_index[perm[a]]
        out[idx] = image_block
    return out


@dataclass(frozen=True)
class Perm:
    """Hashable wrapper for 1-indexed permutation array."""

    p: Tuple[int, ...]

    @property
    def n(self) -> int:
        return len(self.p) - 1


def closure_generated_by(gens: List[List[int]]) -> List[Perm]:
    """Compute group closure by BFS (right-multiplying by generators)."""
    if not gens:
        raise ValueError("Need at least one generator")
    n = len(gens[0]) - 1
    for g in gens:
        if len(g) != n + 1:
            raise ValueError("Generators must have same size")
    gens_t = [Perm(tuple(g)) for g in gens]
    e = Perm(tuple([0] + list(range(1, n + 1))))
    seen = {e.p}
    q = deque([e])
    out = [e]
    while q:
        cur = q.popleft()
        for g in gens_t:
            nxt = compose(list(cur.p), list(g.p))  # cur ∘ g
            t = tuple(nxt)
            if t in seen:
                continue
            seen.add(t)
            p = Perm(t)
            out.append(p)
            q.append(p)
    return out


def extract_tau_sigma(path: Path) -> Optional[Tuple[str, str]]:
    tau = None
    sigma = None
    tau_re = re.compile(r"\btau\s*:=\s*(.*)$")
    sigma_re = re.compile(r"\bsigma\s*:=\s*(.*)$")
    b = path.read_bytes()
    # Many of our PowerShell Tee-Object logs are UTF-16LE with BOM.
    if b.startswith(b"\xff\xfe") or b.startswith(b"\xfe\xff"):
        text = b.decode("utf-16", errors="ignore")
    else:
        text = b.decode("utf-8", errors="ignore")
    for raw in text.splitlines():
        # Remove non-printing / ANSI-ish escapes by keeping basic ASCII + common punctuation
        line = raw.strip()
        m = tau_re.search(line)
        if m:
            tau = m.group(1).strip()
        m = sigma_re.search(line)
        if m:
            sigma = m.group(1).strip()
    if tau and sigma is not None:
        return tau, sigma
    return None


def main() -> None:
    logs_dir = Path(__file__).resolve().parents[1] / "logs"
    candidates = sorted(logs_dir.glob("D24_*.log"))
    if not candidates:
        print(f"No D24 logs found in {logs_dir}")
        return

    n = 24
    print(f"Scanning {len(candidates)} logs in {logs_dir} for tau/sigma...\n")

    for p in candidates:
        ts = extract_tau_sigma(p)
        if not ts:
            continue
        tau_s, sigma_s = ts
        tau = parse_gap_cycles(tau_s, n)
        sigma = parse_gap_cycles(sigma_s, n)
        gens = [tau, sigma]
        group = closure_generated_by(gens)

        orders = [perm_order(g.p) for g in group]
        max_ord = max(orders)
        order_hist = Counter(orders)
        ctau = cycle_type(tau)
        csigma = cycle_type(sigma)
        blocks = block_decomposition_from_tau(tau)
        sigma_blocks = induced_block_permutation(sigma, blocks)
        csigma_blocks = cycle_type(sigma_blocks)

        print(f"== {p.name} ==")
        print(f"- |<tau,sigma>| = {len(group)}")
        print(f"- max element order = {max_ord}")
        print(f"- tau order={perm_order(tau)}, cycle_type={ctau}")
        print(f"- sigma order={perm_order(sigma)}, cycle_type={csigma}")
        print(f"- sigma block-cycle type (S12)={csigma_blocks}")
        print(f"- element-order histogram = {dict(sorted(order_hist.items()))}")
        print()


if __name__ == '__main__':
    main()

