"""Run initial-line saturation and independent verification, preserving timing."""
from pathlib import Path
import datetime
import hashlib
import json
import subprocess
import time

OUT = Path(__file__).resolve().parent
PROJECT = OUT.parents[2]
BOOTSTRAP = Path((OUT.parent / "LATEST_BOOTSTRAP.txt").read_text().strip())
AMBIENT = BOOTSTRAP / "bootstrap" / "ambient.gp"
ORBITS = OUT.parent / "witnesses" / "orbits" / "orbit_representatives.gp"
launcher = '\n'.join([f'AMBIENT={json.dumps(str(AMBIENT))};',
                      f'ORBITS={json.dumps(str(ORBITS))};',
                      f'OUT={json.dumps(str(OUT))};',
                      f'read({json.dumps(str(OUT / "build_initial4.gp"))});', 'quit;'])
(OUT / "launch.gp").write_text(launcher+"\n")
for name in ["FIRST4.gp", "records.jsonl", "ambient.json", "gp_timing.json"]:
    (OUT / name).unlink(missing_ok=True)
timings = {"started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
           "ambient": str(AMBIENT), "orbits": str(ORBITS),
           "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in [AMBIENT, ORBITS]}}
for name, cmd in [("gp", ["/opt/homebrew/bin/gp", "-q", str(OUT / "launch.gp")]),
                  ("independent_verification", ["python3", str(OUT / "verify_initial4.py"), str(OUT)])]:
    start = time.perf_counter()
    p = subprocess.run(cmd, text=True, capture_output=True)
    elapsed = time.perf_counter()-start
    (OUT / f"{name}.stdout").write_text(p.stdout)
    (OUT / f"{name}.stderr").write_text(p.stderr)
    timings[name] = {"wall_seconds": elapsed, "returncode": p.returncode, "command": cmd}
    (OUT / "timings.json").write_text(json.dumps(timings, indent=2)+"\n")
    print(name, p.returncode, elapsed, p.stdout, p.stderr)
    if p.returncode or "***" in p.stderr:
        raise SystemExit("Computation failed; outputs retained.")
timings["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
(OUT / "timings.json").write_text(json.dumps(timings, indent=2)+"\n")
