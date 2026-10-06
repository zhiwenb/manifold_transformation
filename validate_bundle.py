"""Check bundled input integrity without accessing external project files."""
from pathlib import Path
import hashlib, json
root = Path(__file__).resolve().parent
records = json.loads((root / 'provenance/input_hashes.json').read_text())
for item in records:
    path = root / item['file']
    assert path.is_file(), f"Missing input: {item['file']}"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], item['file']
index = json.loads((root / "provenance/input_index.json").read_text())
for item in index:
    assert (root / item["file"]).is_file(), item["file"]
print(f"Verified {len(records)} packaged inputs and {len(index)} plotting views.")
