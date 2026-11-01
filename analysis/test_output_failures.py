"""Exercise output-open and buffered-write failures without touching project data."""
from pathlib import Path
import subprocess
import sys
import tempfile
binary = str(Path(sys.argv[1]).resolve())
with tempfile.TemporaryDirectory(prefix="rutherford-io-") as temp:
    root = Path(temp)
    (root / "data").mkdir()
    target = root / "data/nuclear_timestep_convergence.csv"
    target.mkdir()
    result = subprocess.run([binary, "--convergence"], cwd=root, capture_output=True, text=True)
    assert result.returncode != 0 and "ERROR:" in result.stderr, result
    assert "Saved data/" not in result.stdout
    target.rmdir()
    if Path("/dev/full").exists():
        target.symlink_to("/dev/full")
        result = subprocess.run([binary, "--convergence"], cwd=root, capture_output=True, text=True)
        assert result.returncode != 0 and "ERROR:" in result.stderr, result
        assert "Saved data/" not in result.stdout
print("PASS: output-open and buffered-write failures produce nonzero exit status")
