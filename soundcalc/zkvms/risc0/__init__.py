from soundcalc.zkvms.zkvm import zkVM
from pathlib import Path

def load():
    return zkVM.load_fri_from_toml(Path(__file__).parent / "risc0.toml")
