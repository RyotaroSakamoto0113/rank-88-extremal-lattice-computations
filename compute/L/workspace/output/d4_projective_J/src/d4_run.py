#!/usr/bin/env python3
"""Timed, bounded parallel execution of the 78 rank-four anchor certificates.

Every process is reaped with wait4(its_pid, 0), so CPU times remain accurate
when processes overlap. Only a run covering all 78 anchors can be complete.
Source code and inputs are copied before dispatch; no mutable inputs are read
by the workers. This driver does not execute upon import.
"""
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time

HERE = Path(__file__).resolve().parent
COUNT_NAMES = ["raw_rank2", "nonsaturated_rank2", "ordinary_rank2_rejections",
               "duplicate_rank2", "kept_rank2", "raw_rank3",
               "nonsaturated_rank3", "ordinary_rank3_rejections",
               "second_order_rank3_rejections", "duplicate_rank3",
               "kept_rank3", "raw_rank4", "exact_forms"]


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def copy_snapshot(source, target, executable=False):
    need(source.is_file(), f"missing input: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    target.chmod(0o555 if executable else 0o444)


def prerequisite(pointer, summary_name):
    path = Path((HERE / pointer).read_text().strip()).resolve()
    summary = json.loads((path / summary_name).read_text())
    need(summary.get("status") == "complete", f"incomplete prerequisite: {path / summary_name}")
    return path, summary


def parse_anchors(value):
    if value is None:
        return list(range(78))
    try:
        result = [int(x.strip()) for x in value.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("anchors must be comma-separated integers") from exc
    if not result or len(result) != len(set(result)) or not all(0 <= x < 78 for x in result):
        raise argparse.ArgumentTypeError("anchors must be distinct IDs between 0 and 77")
    return sorted(result)


def pari_errors(path):
    return [line for line in path.read_text(errors="replace").splitlines()
            if line.lstrip().startswith("***") and "Warning:" not in line]


def read_gp_vector(path, name):
    match = re.search(r"\b" + re.escape(name) + r"\s*=\s*(\[[^\[\]]*\])\s*;", path.read_text())
    need(match is not None, f"result.gp has no {name} vector")
    return json.loads(match.group(1))


def read_counts(path):
    values = read_gp_vector(path, "COUNTS")
    need(len(values) >= 13, "COUNTS must contain at least 13 entries")
    need(all(type(n) is int and n >= 0 for n in values[:13]), "invalid COUNTS entry")
    return values


def verify_exact_counts(shard, output_path, expected_forms, counts=None, anchor_result=None):
    with (shard / "exact_manifest.csv").open(newline="") as stream:
        manifest = [row for row in csv.reader(stream) if row]
    exact = [json.loads(line) for line in output_path.read_text().splitlines() if line.strip()]
    need(len(exact) == len(manifest) == expected_forms > 0, "GP/manifest/exact batch coverage mismatch")
    signed_vectors = 0
    nodes = 0
    ldl_seconds = 0.0
    enumeration_seconds = 0.0
    by_dimension = {}
    for number, (row, record) in enumerate(zip(manifest, exact), 1):
        need(len(row) >= 4, f"short manifest row {number}")
        fid, dimension, bound, count = map(int, row[:4])
        need(fid == number, f"manifest form number mismatch at {number}")
        need(int(record["file"]) == number, f"exact form number mismatch at {number}")
        need(record["dimension"] == dimension and int(record["bound"]) == bound,
             f"exact dimension/bound mismatch at form {number}")
        need(record["complete"] is True, f"incomplete exact enumeration at form {number}")
        need(record["signed_vectors"] == count and 2 * record["pairs"] == count,
             f"signed vector count mismatch at form {number}")
        signed_vectors += count
        nodes += record.get("nodes", 0)
        ldl_seconds += record.get("ldl_seconds", 0)
        enumeration_seconds += record.get("enumeration_seconds", 0)
        stats = by_dimension.setdefault(dimension, {"forms": 0, "signed_vectors": 0})
        stats["forms"] += 1
        stats["signed_vectors"] += count
        if counts is not None:
            need(dimension in (20, 18, 16), f"unexpected quotient dimension {dimension}")
            need(len(row) >= 5 and row[4] == {20: "second", 18: "third", 16: "fourth"}[dimension],
                 f"stage/dimension mismatch at form {number}")
    if counts is not None:
        need(anchor_result is not None and len(anchor_result) >= 2, "missing anchor norm")
        m = anchor_result[1]
        need(counts[0] == sum(counts[1:5]), "rank-two disposition counts do not partition raw candidates")
        need(counts[5] == sum(counts[6:11]), "rank-three disposition counts do not partition raw candidates")
        expected_dimension_forms = {20: 3, 18: 3 * counts[4], 16: 3 * counts[10]}
        expected_dimension_signed = {20: 2 * counts[0], 18: 2 * counts[5], 16: 2 * counts[11]}
        for dimension in [20, 18, 16]:
            stats = by_dimension.get(dimension, {"forms": 0, "signed_vectors": 0})
            need(stats["forms"] == expected_dimension_forms[dimension],
                 f"{dimension}-dimensional form coverage mismatch")
            need(stats["signed_vectors"] == expected_dimension_signed[dimension],
                 f"{dimension}-dimensional signed count differs from raw candidate counter")
        need(sum(expected_dimension_forms.values()) == counts[12], "form total differs from COUNTS")
    return {"forms": len(exact), "signed_vectors": signed_vectors, "nodes": nodes,
            "sum_ldl_seconds": ldl_seconds, "sum_enumeration_seconds": enumeration_seconds,
            "all_manifest_counts_match": True, "by_dimension": by_dimension,
            "all_stage_partitions_match": counts is not None}


class Runner:
    def __init__(self, root, gp_command):
        self.root = root
        self.inputs = root / "inputs"
        self.gp_command = gp_command
        self.children = {}
        self.lock = threading.Lock()
        self.stop = threading.Event()

    def terminate_children(self):
        self.stop.set()
        with self.lock:
            for pid in list(self.children):
                try:
                    os.kill(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass

    def process(self, shard, name, command):
        started = utcnow()
        wall = time.perf_counter()
        command = list(map(str, command))
        stage = {"stage": name, "command": command, "started_utc": started,
                 "stdout": str(shard / (name + ".stdout")),
                 "stderr": str(shard / (name + ".stderr"))}
        save_json(shard / (name + ".command.json"), stage)
        with (shard / (name + ".stdout")).open("w") as stdout, \
             (shard / (name + ".stderr")).open("w") as stderr:
            with self.lock:
                need(not self.stop.is_set(), "run cancelled before process start")
                process = subprocess.Popen(command, cwd=shard, stdout=stdout, stderr=stderr)
                self.children[process.pid] = name
                stage["pid"] = process.pid
                save_json(shard / (name + ".command.json"), stage)
                save_json(shard / (name + ".process.json"),
                          {"pid": process.pid, "stage": name, "status": "running",
                           "started_utc": started, "command": command})
            try:
                while True:
                    try:
                        pid, status, usage = os.wait4(process.pid, 0)
                        break
                    except InterruptedError:
                        continue
                need(pid == process.pid, "wait4 returned wrong process")
                process.returncode = os.waitstatus_to_exitcode(status)
            finally:
                with self.lock:
                    self.children.pop(process.pid, None)
        stage.update(returncode=process.returncode, finished_utc=utcnow(),
                     wall_seconds=time.perf_counter() - wall,
                     user_cpu_seconds=usage.ru_utime, system_cpu_seconds=usage.ru_stime,
                     max_rss_native_units=usage.ru_maxrss,
                     voluntary_context_switches=usage.ru_nvcsw,
                     involuntary_context_switches=usage.ru_nivcsw)
        save_json(shard / (name + ".timing.json"), stage)
        save_json(shard / (name + ".process.json"),
                  {"pid": stage["pid"], "stage": name, "status": "finished",
                   "started_utc": started, "finished_utc": stage["finished_utc"],
                   "returncode": stage["returncode"], "command": command})
        return stage

    def anchor(self, anchor):
        shard = self.root / "anchors" / f"{anchor:03d}"
        shard.mkdir(parents=True)
        record = {"anchor_id": anchor, "status": "running", "started_utc": utcnow(),
                  "directory": str(shard), "stages": []}
        started = time.perf_counter()
        try:
            for name in ["ambient.gp", "generators.gp"]:
                copy_snapshot(self.inputs / name, shard / name)
            launcher = shard / "launch.gp"
            lines = ["default(parisizemax,3000000000);",
                     "default(nbthreads,1);",
                     f"OUT={json.dumps(str(shard))};ANCHOR={anchor};SUPPORT={json.dumps(str(self.inputs))};"]
            for path in [self.root / "src" / "common.gp",
                         self.inputs / "ordinary_prune" / "ordinary_lookup.gp",
                         self.inputs / "ordinary_prune" / "shorter_ordinary.gp",
                         self.inputs / "FIRST4.gp", self.inputs / "STAB4.gp",
                         self.inputs / "ordinary_representatives.gp",
                         self.root / "src" / "d4.gp"]:
                lines.append(f"read({json.dumps(str(path))});")
            lines.append("quit;")
            launcher.write_text("\n".join(lines) + "\n")
            launcher.chmod(0o444)
            record["input_hashes"] = {name: digest(shard / name)
                                      for name in ["ambient.gp", "generators.gp", "launch.gp"]}
            record["shared_input_manifest"] = str(self.root / "INPUTS_SHA256.json")
            save_json(shard / "ANCHOR_SUMMARY.json", record)
            stage = self.process(shard, "gp_generation", [self.gp_command, "-fq", launcher])
            record["stages"].append(stage)
            save_json(shard / "ANCHOR_SUMMARY.json", record)
            need(stage["returncode"] == 0, "GP process failed")
            need(not pari_errors(shard / "gp_generation.stderr"), "GP reported an error")
            need("D4_ANCHOR_COMPLETE" in (shard / "gp_generation.stdout").read_text(),
                 "GP anchor completion marker missing")
            for name in ["exact_inputs.txt", "exact_manifest.csv", "enumerated_vectors.gp", "result.gp"]:
                need((shard / name).is_file(), f"GP omitted {name}")
            values = read_counts(shard / "result.gp")
            anchor_result = read_gp_vector(shard / "result.gp", "ANCHOR_RESULT")
            need(len(anchor_result) >= 5 and all(type(v) is int for v in anchor_result[:5]),
                 "invalid ANCHOR_RESULT vector")
            need(anchor_result[0] == anchor, "result belongs to a different anchor")
            reps_text = (self.inputs / "ordinary_representatives.gp").read_text().strip()
            reps = json.loads(reps_text.removeprefix("MINIMA_ORBIT_REPS=").removesuffix(";"))
            need(reps[anchor][0] == anchor and anchor_result[1] == reps[anchor][1],
                 "anchor norm differs from original ordinary orbit representative")
            record["anchor_result"] = anchor_result
            record["counts"] = dict(zip(COUNT_NAMES, values[:13]))
            record["raw_counts_vector"] = values
            record["potential_counterexample"] = values[11] > 0
            save_json(shard / "ANCHOR_SUMMARY.json", record)
            stage = self.process(shard, "exact_enumeration",
                                 [self.inputs / "exact_enum", "--batch", shard / "exact_inputs.txt", "60"])
            record["stages"].append(stage)
            save_json(shard / "ANCHOR_SUMMARY.json", record)
            need(stage["returncode"] == 0, "independent exact enumeration process failed")
            record["independent_verification"] = verify_exact_counts(
                shard, shard / "exact_enumeration.stdout", values[12], values, anchor_result)
            record["status"] = "counterexample" if values[11] else "complete"
            record["output_hashes"] = {name: digest(shard / name) for name in
                                       ["exact_inputs.txt", "exact_manifest.csv", "enumerated_vectors.gp",
                                        "result.gp", "exact_enumeration.stdout"]}
        except BaseException as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        record.update(finished_utc=utcnow(), wall_seconds=time.perf_counter() - started,
                      user_cpu_seconds=sum(s["user_cpu_seconds"] for s in record["stages"]),
                      system_cpu_seconds=sum(s["system_cpu_seconds"] for s in record["stages"]))
        save_json(shard / "ANCHOR_SUMMARY.json", record)
        return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", type=parse_anchors, help="comma-separated subset, default all IDs 0..77")
    parser.add_argument("--jobs", type=int, default=4, help="maximum concurrent anchor processes, default 4")
    parser.add_argument("--label", default="", help="optional run-directory label")
    parser.add_argument("--gp", default="gp", help="GP executable, default gp")
    args = parser.parse_args()
    anchors = args.anchors if args.anchors is not None else list(range(78))
    need(args.jobs > 0, "--jobs must be positive")
    need(hasattr(os, "wait4"), "per-process CPU measurement requires os.wait4")
    gp_command = shutil.which(args.gp)
    need(gp_command is not None, "GP executable was not found")
    boot, boot_summary = prerequisite("LATEST_BOOTSTRAP.txt", "BOOTSTRAP_SUMMARY.json")
    d3, d3_summary = prerequisite("LATEST_D3.txt", "D3_SUMMARY.json")
    ordinary, ordinary_summary = prerequisite("LATEST_ORDINARY.txt", "ORDINARY_SUMMARY.json")
    need(Path(d3_summary["bootstrap_directory"]).resolve() == boot, "d3/bootstrap input mismatch")
    need(ordinary_summary.get("ordinary_orbits") == 78, "ordinary prerequisite does not cover 78 orbits")
    need(ordinary_summary.get("previous_lookup_data_identical") is True,
         "ordinary lookup has not been matched to the fresh complete sphere")
    need((HERE / "src" / "d4.gp").is_file(), "src/d4.gp is not ready")
    label = re.sub(r"[^A-Za-z0-9_-]", "_", args.label)
    prefix = datetime.datetime.now().strftime("%Y%m%d-%H%M%S-d4-") + (label + "-" if label else "")
    (HERE / "runs").mkdir(exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix=prefix, dir=HERE / "runs"))
    for path in (HERE / "src").glob("*.gp"):
        copy_snapshot(path, root / "src" / path.name)
    copy_snapshot(Path(__file__), root / "src" / "d4_run.py")
    for name in ["ambient.gp", "generators.gp"]:
        copy_snapshot(boot / "bootstrap" / name, root / "inputs" / name)
    for name in ["FIRST4.gp", "STAB4.gp"]:
        copy_snapshot(HERE / "initial4" / name, root / "inputs" / name)
    copy_snapshot(HERE / "witnesses" / "orbits" / "orbit_representatives.gp",
                  root / "inputs" / "ordinary_representatives.gp")
    ordinary_prune = HERE / "ordinary_prune"
    if not ordinary_prune.is_dir():
        ordinary_prune = HERE / "witnesses" / "ordinary_prune"
    for name in ["ordinary_lookup.gp", "shorter_ordinary.gp"]:
        copy_snapshot(ordinary_prune / name, root / "inputs" / "ordinary_prune" / name)
    copy_snapshot(boot / "exact_enum", root / "inputs" / "exact_enum", executable=True)
    if (boot / "src" / "exact_enum.cpp").is_file():
        copy_snapshot(boot / "src" / "exact_enum.cpp", root / "src" / "exact_enum.cpp")
    for name, value in [("bootstrap", boot_summary), ("d3", d3_summary), ("ordinary", ordinary_summary)]:
        save_json(root / (name + "_prerequisite.json"), value)
    hashes = {str(path.relative_to(root)): digest(path)
              for directory in [root / "src", root / "inputs"]
              for path in sorted(directory.rglob("*")) if path.is_file()}
    save_json(root / "INPUTS_SHA256.json", hashes)
    started = time.perf_counter()
    summary = {"status": "running", "run_directory": str(root), "started_utc": utcnow(),
               "requested_anchors": anchors, "jobs": args.jobs, "platform": sys.platform,
               "bootstrap_directory": str(boot), "d3_directory": str(d3),
               "ordinary_directory": str(ordinary), "input_hashes": hashes, "anchors": {}}
    runner = Runner(root, gp_command)
    pool = ThreadPoolExecutor(max_workers=args.jobs)
    pending = {}
    remaining = iter(anchors)
    had_error = False
    interrupted = None

    def save():
        verified = [r for r in summary["anchors"].values() if r["status"] in ("complete", "counterexample")]
        summary["counts"] = {name: sum(r["counts"][name] for r in verified) for name in COUNT_NAMES}
        summary["completed_anchor_ids"] = sorted(r["anchor_id"] for r in verified)
        summary["wall_seconds"] = time.perf_counter() - started
        summary["user_cpu_seconds"] = sum(r["user_cpu_seconds"] for r in summary["anchors"].values())
        summary["system_cpu_seconds"] = sum(r["system_cpu_seconds"] for r in summary["anchors"].values())
        save_json(root / "D4_SUMMARY.json", summary)

    def submit_one():
        anchor = next(remaining, None)
        if anchor is not None:
            pending[pool.submit(runner.anchor, anchor)] = anchor

    print("RUN_DIRECTORY", root, flush=True)
    (HERE / "LATEST_D4_RUN.txt").write_text(str(root) + "\n")
    save()
    try:
        for _ in range(min(args.jobs, len(anchors))):
            submit_one()
        while pending:
            finished, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in finished:
                anchor = pending.pop(future)
                result = future.result()
                summary["anchors"][str(anchor)] = result
                if result["status"] != "complete":
                    had_error = True
                print("ANCHOR", anchor, result["status"],
                      f"wall={result['wall_seconds']:.6f}",
                      "counts=" + json.dumps(result.get("counts", {}), separators=(",", ":")), flush=True)
            save()
            if not had_error:
                for _ in finished:
                    submit_one()
    except BaseException as exc:
        interrupted = f"{type(exc).__name__}: {exc}"
        runner.terminate_children()
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
        for future, anchor in pending.items():
            if future.done() and not future.cancelled():
                try:
                    summary["anchors"][str(anchor)] = future.result()
                except BaseException as exc:
                    interrupted = interrupted or f"{type(exc).__name__}: {exc}"
    all_complete = (set(summary["anchors"]) == {str(n) for n in range(78)}
                    and all(r["status"] == "complete" for r in summary["anchors"].values()))
    counterexample = any(r["status"] == "counterexample" for r in summary["anchors"].values())
    failed = interrupted or any(r["status"] == "failed" for r in summary["anchors"].values())
    summary.update(status="failed" if failed else "counterexample" if counterexample else "complete" if all_complete else "partial",
                   finished_utc=utcnow(), exclusion_complete=all_complete,
                   verified_counterexample=counterexample,
                   not_started_anchor_ids=sorted(set(anchors) - {int(n) for n in summary["anchors"]}))
    if interrupted:
        summary["error"] = interrupted
    save()
    if all_complete:
        (HERE / "LATEST_D4.txt").write_text(str(root) + "\n")
    print("D4_RUN_FINISHED", summary["status"], "wall", summary["wall_seconds"],
          "user_cpu", summary["user_cpu_seconds"], "system_cpu", summary["system_cpu_seconds"], flush=True)
    if summary["status"] not in ("complete", "partial"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
