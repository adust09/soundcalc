from __future__ import annotations

from soundcalc.zkvms.fri_based_vm import FRIBasedCircuit, FRIBasedVM, FRIBasedVMConfig

from ..common.fields import *


def _parse_fri_folding_arities(fri_fa_str: str) -> tuple[list[int], int]:
    """
    Parse FRI folding arity string like "23-19-15-11-8-5" into folding factors and early stop degree.

    The sequence represents log2 of domain sizes at each step.
    E.g., "23-19-15-11-8-5" means:
      - 23→19: fold by 2^(23-19) = 16
      - 19→15: fold by 2^4 = 16
      - 15→11: fold by 2^4 = 16
      - 11→8: fold by 2^3 = 8
      - 8→5: fold by 2^3 = 8
      - Final 5: early_stop_degree = 2^5 = 32
    """
    values = [int(x) for x in fri_fa_str.split('-')]
    folding_factors = []
    for i in range(len(values) - 1):
        fold_log = values[i] - values[i + 1]
        folding_factors.append(1 << fold_log)
    early_stop_degree = 1 << values[-1]
    return folding_factors, early_stop_degree


def _make_circuit(name, bits, bf, d, fixed, stage1, pols, queries, opens, fri_fa) -> FRIBasedCircuit:
    """Factory function to create a circuit from table parameters."""
    FRI_folding_factors, FRI_early_stop_degree = _parse_fri_folding_arities(fri_fa)

    return FRIBasedCircuit(FRIBasedVMConfig(
        name=name,
        trace_length=1 << bits,
        rho=1 / (1 << bf),
        AIR_max_degree=d,
        num_columns=fixed + stage1,
        batch_size=pols,
        num_queries=queries,
        max_combo=opens,
        FRI_folding_factors=FRI_folding_factors,
        FRI_early_stop_degree=FRI_early_stop_degree,
        field=GOLDILOCKS_3,
        hash_size_bits=256,
        power_batching=True,
        grinding_query_phase=0,
    ))


# Base ZisK circuits
# Columns: (name, bits, bf, d, fixed, stage1, pols, queries, opens, fri_fa)
#
# Aligned for better reading:
#                              bits bf  d  fix stg1  pols  qry opn  fri_fa
ZISK_BASE_CIRCUITS = [
    ("Main",                    22,  1, 3,   3,  38,   61, 128,  3, "23-19-15-11-8-5"),
    ("Rom",                     22,  1, 2,   1,   1,   18, 128,  3, "23-19-15-11-8-5"),
]


class ZiskPreset:
    @staticmethod
    def default() -> FRIBasedVM:
        """
        Create a ZisK VM with multiple circuits.

        For ZisK, we populate the trace parameters from its constraint description:
            https://github.com/0xPolygonHermez/zisk/blob/main/pil/zisk.pil

        The rest of the parameters are adapted from the "eSTARK: Extending STARKs with Arguments" paper:
           https://eprint.iacr.org/2023/474
        """
        circuits = [_make_circuit(*row) for row in ZISK_BASE_CIRCUITS]
        return FRIBasedVM(name="ZisK", circuits=circuits)


if __name__ == "__main__":
    import unittest

    class TestParseFriFoldingArities(unittest.TestCase):
        def test_main_circuit_fri_fa(self):
            # 23-19-15-11-8-5: folds by 16,16,16,8,8 then stops at degree 32
            factors, early_stop = _parse_fri_folding_arities("23-19-15-11-8-5")
            self.assertEqual(factors, [16, 16, 16, 8, 8])
            self.assertEqual(early_stop, 32)

        def test_short_fri_fa(self):
            # 17-13-9-5: folds by 16,16,16 then stops at degree 32
            factors, early_stop = _parse_fri_folding_arities("17-13-9-5")
            self.assertEqual(factors, [16, 16, 16])
            self.assertEqual(early_stop, 32)

        def test_final_circuit_fri_fa(self):
            # 20-15-10: folds by 32,32 then stops at degree 1024
            factors, early_stop = _parse_fri_folding_arities("20-15-10")
            self.assertEqual(factors, [32, 32])
            self.assertEqual(early_stop, 1024)

    unittest.main()
