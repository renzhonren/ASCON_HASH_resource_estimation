"""
The function layout follows ASCON_quantum-master/.../ASCON_HASH.py: each round
is built from constant addition, a substitution layer, and a linear layer.
Use --scheme to select QAND out-of-place, Amy out-of-place, or Amy in-place.
"""

import argparse
import json
from functools import lru_cache
from pathlib import Path
import re

from projectq import MainEngine
from projectq.backends import ResourceCounter
from projectq.meta import Dagger
from projectq.ops import CNOT, H, Measure, S, T, Tdag, X


ROUND_CONSTANTS = (0xF0, 0xE1, 0xD2, 0xC3, 0xB4, 0xA5,
                   0x96, 0x87, 0x78, 0x69, 0x5A, 0x4B)
ABSORB_ROUNDS = (1, 13, 25, 37)
PADDING_ROUND = 49  # four data blocks, then a separate 0x80 block
SQUEEZE_OUTPUT_ROUNDS = (60, 72, 84, 96)
PRECOMPUTED_INITIAL_STATE = (
    0xEE9398AADB67F03D, 0x8BB21831C60F1002,
    0xB48A92DB98D5DA62, 0x43189921B8F8E3E8,
    0x348FA5C9D525E140,
)


@lru_cache(maxsize=1)
def Linear_CNOT_data():
    """Read Roy et al.'s 1595-CNOT, depth-119 in-place linear circuit."""
    path = Path(__file__).with_name("roy_linear_xzlbz.txt")
    if not path.is_file():
        raise FileNotFoundError(
            "Download Roy's implementations/xzlbz.txt as "
            "roy_linear_xzlbz.txt in this directory; see README.md")
    lines = path.read_text(encoding="utf-8").splitlines()
    permutation = [int(v) for v in re.findall(r"x\[(\d+)\]", lines[0])]
    pattern = re.compile(r"x\[(\d+)\] = x\[\d+\] \^ x\[(\d+)\]")
    cnots = [tuple(map(int, pattern.fullmatch(line).groups())) for line in lines[1:]]
    if sorted(permutation) != list(range(320)) or len(cnots) != 1595:
        raise ValueError("Roy linear-layer data are incomplete")
    return permutation, cnots


def Toffoli_gate(eng, a, b, c, method="amy", work=None):
    """Amy Figure 2 decomposition or clean-target QAND from Figure 3."""
    if method == "qand":
        if work is None:
            raise ValueError("QAND requires a clean work qubit")
        H | c
        CNOT | (b, work)
        CNOT | (c, a)
        CNOT | (c, b)
        CNOT | (a, work)
        Tdag | a
        Tdag | b
        T | c
        T | work
        CNOT | (a, work)
        CNOT | (c, b)
        CNOT | (c, a)
        CNOT | (b, work)
        H | c
        S | c
    elif method == "amy":
        # Same gate order as the reference ASCON_HASH.py, verified as Toffoli.
        Tdag | a
        Tdag | b
        H | c
        CNOT | (c, a)
        T | a
        CNOT | (b, c)
        CNOT | (b, a)
        T | c
        Tdag | a
        CNOT | (b, c)
        CNOT | (c, a)
        T | a
        Tdag | c
        CNOT | (b, a)
        H | c
    else:
        raise ValueError(method)


def QAND_dagger_cost(eng, a, b, c):
    """Charge the measured correction branch used in the paper's QAND model."""
    H | c
    H | b
    Measure | c
    CNOT | (a, b)
    X | c
    H | b


def Message_XOR(eng, message, x0, block):
    for i in range(64):
        CNOT | (message[64 * block + i], x0[i])


def S_constant_XOR(eng, x2, value):
    # Roy's linear network uses MSB-first wires; the 8-bit constant is low.
    for bit in range(8):
        if value & (1 << (7 - bit)):
            X | x2[56 + bit]


def add_constant(eng, x2, round_number):
    S_constant_XOR(eng, x2, ROUND_CONSTANTS[(round_number - 1) % 12])


def Substitution_Layer(eng, x0, x1, x2, x3, x4,
                       top0, top1, top2, top3, top4,
                       new0, new1, new2, new3, new4,
                       method="amy", work=None):
    """Figure 7 out-of-place S-boxes, 64 in parallel; top bits start at one."""
    x = (x0, x1, x2, x3, x4)
    top = (top0, top1, top2, top3, top4)
    new = (new0, new1, new2, new3, new4)
    for i in range(64):
        CNOT | (x4[i], x0[i])
        CNOT | (x1[i], x2[i])
        CNOT | (x3[i], x4[i])
        for word in range(5):
            CNOT | (x[word][i], top[word][i])
        for output, left, right in ((0, 1, 2), (1, 2, 3), (2, 3, 4),
                                    (3, 4, 0), (4, 0, 1)):
            ancilla = work[output][i] if work is not None else None
            Toffoli_gate(eng, top[left][i], x[right][i], new[output][i],
                         method, ancilla)
        for word in range(5):
            CNOT | (x[word][i], new[word][i])
        CNOT | (new0[i], new1[i])
        CNOT | (new2[i], new3[i])
        CNOT | (new4[i], new0[i])
        X | new2[i]
        for word in range(5):
            CNOT | (x[word][i], top[word][i])


def Substitution_Layer_dag(eng, x, top, new, method="amy", work=None):
    if method == "amy":
        with Dagger(eng):
            Substitution_Layer(eng, *x, *top, *new, method="amy")
        return

    # QAND-dagger uses a measured correction and cannot use Dagger(eng).
    for i in range(64):
        for word in reversed(range(5)):
            CNOT | (x[word][i], top[word][i])
        X | new[2][i]
        CNOT | (new[4][i], new[0][i])
        CNOT | (new[2][i], new[3][i])
        CNOT | (new[0][i], new[1][i])
        for word in reversed(range(5)):
            CNOT | (x[word][i], new[word][i])
        for output, left, right in ((4, 0, 1), (3, 4, 0), (2, 3, 4),
                                    (1, 2, 3), (0, 1, 2)):
            QAND_dagger_cost(eng, top[left][i], x[right][i], new[output][i])
        for word in reversed(range(5)):
            CNOT | (x[word][i], top[word][i])
        CNOT | (x[3][i], x[4][i])
        CNOT | (x[1][i], x[2][i])
        CNOT | (x[4][i], x[0][i])


def Inplace_Substitution_Layer(eng, x0, x1, x2, x3, x4):
    """Figure 6 five-wire S-box, directly on the five 64-bit state words."""
    for i in range(64):
        CNOT | (x4[i], x0[i])
        CNOT | (x1[i], x2[i])
        CNOT | (x3[i], x4[i])
        X | x1[i]
        X | x3[i]
        Toffoli_gate(eng, x1[i], x2[i], x0[i])
        Toffoli_gate(eng, x3[i], x4[i], x2[i])
        Toffoli_gate(eng, x1[i], x2[i], x0[i])
        for bit in (x0[i], x1[i], x2[i], x3[i]):
            X | bit
        Toffoli_gate(eng, x0[i], x1[i], x4[i])
        X | x0[i]
        X | x4[i]
        Toffoli_gate(eng, x2[i], x3[i], x1[i])
        X | x1[i]
        X | x2[i]
        Toffoli_gate(eng, x0[i], x4[i], x3[i])
        Toffoli_gate(eng, x1[i], x2[i], x0[i])
        X | x1[i]
        X | x4[i]
        CNOT | (x0[i], x1[i])
        CNOT | (x2[i], x3[i])
        CNOT | (x4[i], x0[i])
        X | x2[i]


def LinearDiffusion_Layer(eng, state, inverse=False):
    permutation, cnots = Linear_CNOT_data()
    wires = [qubit for word in state for qubit in word]
    if inverse:
        for target, control in reversed(cnots):
            CNOT | (wires[control], wires[target])
        wires = [wires[permutation[i]] for i in range(320)]
    else:
        inverse_permutation = [0] * 320
        for index, value in enumerate(permutation):
            inverse_permutation[value] = index
        wires = [wires[inverse_permutation[i]] for i in range(320)]
        for target, control in cnots:
            CNOT | (wires[control], wires[target])
    return [wires[64 * word:64 * (word + 1)] for word in range(5)]


def New_State(eng):
    return [list(eng.allocate_qureg(64)) for _ in range(5)]


def Initialize_Precomputed_State(eng, state):
    """Prepare the fixed state used by the reference HASH implementation."""
    for word, value in enumerate(PRECOMPUTED_INITIAL_STATE):
        for bit in range(64):
            if value & (1 << (63 - bit)):
                X | state[word][bit]


def Permutation_a(eng, round_number, old, new, top, method, work):
    add_constant(eng, old[2], round_number)
    Substitution_Layer(eng, *old, *top, *new, method=method, work=work)
    return LinearDiffusion_Layer(eng, new)


def Permutation_a_dag(eng, round_number, old, new, top, method, work):
    new = LinearDiffusion_Layer(eng, new, inverse=True)
    Substitution_Layer_dag(eng, old, top, new, method=method, work=work)
    add_constant(eng, old[2], round_number)
    return new


def Partial_Recovery(eng, old, round_number):
    for i in range(64):
        CNOT | (old[3][i], old[4][i])
        CNOT | (old[1][i], old[2][i])
        CNOT | (old[4][i], old[0][i])
    add_constant(eng, old[2], round_number)


def Resource_Summary(resource, scheme):
    count = {"t_gates": 0, "cnot_gates": 0,
             "one_qubit_clifford": 0, "measurements": 0}
    for (gate, controls), amount in resource.gate_counts.items():
        name = type(gate).__name__
        if gate == T or gate == Tdag:
            count["t_gates"] += amount
        elif name == "XGate" and controls == 1:
            count["cnot_gates"] += amount
        elif name in ("XGate", "HGate", "SGate", "SdagGate") and controls == 0:
            count["one_qubit_clifford"] += amount
        elif name == "MeasureGate":
            count["measurements"] += amount
    return {"scheme": scheme,
            "qubits": resource.max_width - 512,
            "total_depth": resource.depth_of_dag,
            "physical_peak_qubits_including_message_and_digest": resource.max_width,
            **count}


def Estimate_Inplace(include_initialization=False):
    resource = ResourceCounter()
    eng = MainEngine(resource, engine_list=[])
    message = list(eng.allocate_qureg(256))
    digest = list(eng.allocate_qureg(256))
    state = New_State(eng)
    if include_initialization:
        Initialize_Precomputed_State(eng, state)
    for r in range(1, 97):
        if r in ABSORB_ROUNDS:
            Message_XOR(eng, message, state[0], ABSORB_ROUNDS.index(r))
        if r == PADDING_ROUND:
            X | state[0][0]  # MSB-first representation of 0x8000...0000
        add_constant(eng, state[2], r)
        Inplace_Substitution_Layer(eng, *state)
        state = LinearDiffusion_Layer(eng, state)
        if r in SQUEEZE_OUTPUT_ROUNDS:
            block = SQUEEZE_OUTPUT_ROUNDS.index(r)
            for bit in range(64):
                CNOT | (state[0][bit], digest[64 * block + bit])
    eng.flush()
    result = Resource_Summary(resource, "amy_inplace")
    if include_initialization:
        result["fixed_initial_state_preparation_included"] = True
    return result


def Estimate_Outplace(method="amy", shallow=False,
                      include_initialization=False):
    resource = ResourceCounter()
    eng = MainEngine(resource, engine_list=[])
    message = list(eng.allocate_qureg(256))
    digest = list(eng.allocate_qureg(256))
    states = {0: New_State(eng)}
    if include_initialization:
        Initialize_Precomputed_State(eng, states[0])
    free = []
    top = [New_State(eng) for _ in range(2 if shallow else 1)]
    work = New_State(eng) if method == "qand" else None
    for line in top:
        for word in line:
            for qubit in word:
                X | qubit

    def forward(r, line):
        old = states[r - 1]
        if r in ABSORB_ROUNDS:
            Message_XOR(eng, message, old[0], ABSORB_ROUNDS.index(r))
        if r == PADDING_ROUND:
            X | old[0][0]
        new = free.pop() if free else New_State(eng)
        states[r] = Permutation_a(eng, r, old, new, top[line], method, work)
        if r in SQUEEZE_OUTPUT_ROUNDS:
            block = SQUEEZE_OUTPUT_ROUNDS.index(r)
            for bit in range(64):
                CNOT | (states[r][0][bit], digest[64 * block + bit])

    def inverse(r, line):
        old = states[r - 1]
        new = Permutation_a_dag(eng, r, old, states[r], top[line], method, work)
        if r in (13, 25):
            Message_XOR(eng, message, old[0], ABSORB_ROUNDS.index(r))
        if r == PADDING_ROUND:
            X | old[0][0]
        free.append(new)
        del states[r]

    if shallow:
        for r in range(1, 10):
            forward(r, 0)
        Partial_Recovery(eng, states[8], 9)
        for group in range(2, 12):
            start = 9 * (group - 1) + 1
            end = min(9 * group, 96)
            forwards = list(range(start, end + 1))
            inverses = list(range(start - 2, start - 10, -1))
            if group == 11:
                inverses = inverses[:len(forwards) - 1]
            for step in range(max(len(forwards), len(inverses))):
                if step < len(forwards):
                    forward(forwards[step], (group - 1) % 2)
                if step < len(inverses):
                    inverse(inverses[step], (group - 2) % 2)
            if end % 9 == 0:
                Partial_Recovery(eng, states[end - 1], end)
    else:
        for start in range(1, 91, 9):
            end = start + 8
            for r in range(start, end + 1):
                forward(r, 0)
            Partial_Recovery(eng, states[end - 1], end)
            inverse_rounds = list(range(end - 1, start - 1, -1))
            if start == 82:
                # Six free registers are enough for final rounds 91..96.
                inverse_rounds = inverse_rounds[:6]
            for r in inverse_rounds:
                inverse(r, 0)
        for r in range(91, 97):
            forward(r, 0)

    eng.flush()
    label = "qand_outplace_standard" if method == "qand" else "amy_outplace_shallow"
    result = Resource_Summary(resource, label)
    result["full_inverse_rounds"] = 77 if shallow else 78
    if include_initialization:
        result["fixed_initial_state_preparation_included"] = True
    return result


def main():
    parser = argparse.ArgumentParser(description="ProjectQ ASCON-HASH resource estimates")
    parser.add_argument("--scheme", choices=("all", "qand", "amy-outplace", "amy-inplace"),
                        default="all")
    parser.add_argument("--include-fixed-initial-state", action="store_true",
                        help="Count 150 X gates preparing the precomputed initial state")
    args = parser.parse_args()
    scheme = args.scheme
    result = {}
    if scheme in ("all", "qand"):
        result["qand"] = Estimate_Outplace(method="qand",
                                           include_initialization=args.include_fixed_initial_state)
    if scheme in ("all", "amy-outplace"):
        result["amy-outplace"] = Estimate_Outplace(
            method="amy", shallow=True,
            include_initialization=args.include_fixed_initial_state)
    if scheme in ("all", "amy-inplace"):
        result["amy-inplace"] = Estimate_Inplace(
            include_initialization=args.include_fixed_initial_state)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
