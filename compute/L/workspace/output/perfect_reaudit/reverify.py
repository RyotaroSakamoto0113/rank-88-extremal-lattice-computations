#!/usr/bin/env python3
"""Reverify the saved field input and perfection certificate using only Python.

All fresh results are written to a new directory. The original input, run,
certificate and verification records are never rewritten. Python >= 3.9.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable_recorded_path(text):
    """Map an archived absolute audit path into this extracted audit tree."""
    p = Path(text.strip())
    if not p.is_absolute():
        return (HERE / p).resolve()
    marker = '/output/perfect_reaudit/'
    if marker in str(p):
        return (HERE / str(p).split(marker, 1)[1]).resolve()
    return p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, help='Saved certificate.json; default uses LATEST_CERTIFICATE.txt')
    parser.add_argument('--output', type=Path, help='New results directory; default creates a fresh run')
    args = parser.parse_args()
    cert = args.certificate.resolve() if args.certificate else portable_recorded_path((HERE / 'LATEST_CERTIFICATE.txt').read_text()) / 'certificate.json'
    need(cert.is_file(), 'certificate not found: ' + str(cert))
    need((cert.parent / 'lattice.json').is_file(), 'certificate run lacks its frozen lattice.json')
    need(sha(HERE / 'input/lattice.json') == sha(cert.parent / 'lattice.json'), 'field input differs from the certificate run input')
    if args.output:
        out = args.output.resolve()
        out.mkdir(parents=True, exist_ok=False)
    else:
        (HERE / 'runs').mkdir(exist_ok=True)
        out = Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-reverify-'), dir=HERE / 'runs'))
    started = time.perf_counter()
    report = {'status': 'running', 'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'output': str(out), 'source_certificate': str(cert), 'stages': [],
              'external_prerequisite': 'The previously proved minimum of this lattice is 8.',
              'scope': 'Independent field/input reconstruction and independent split-certificate verification; no LLL search.'}

    def save():
        (out / 'REVERIFICATION.json').write_text(json.dumps(report, indent=2) + '\n')

    print('RUN_DIRECTORY', out, flush=True)
    save()
    try:
        field = out / 'input'
        field.mkdir()
        sources = [HERE / 'input/verify_tensor_input.py', HERE / 'input/lattice.json',
                   *sorted((HERE / 'input').glob('*.txt'))]
        report['source_hashes'] = {str(p): sha(p) for p in sources}
        for path in sources:
            shutil.copy2(path, field / path.name)
        shutil.copy2(cert, out / 'certificate.json')
        verify_source = HERE / 'src/verify_split_certificate.py'
        shutil.copy2(verify_source, out / verify_source.name)
        report['source_hashes'][str(cert)] = sha(cert)
        report['source_hashes'][str(verify_source)] = sha(verify_source)

        def run(name, command):
            t = time.perf_counter()
            record = {'name': name, 'command': [str(v) for v in command],
                      'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
            report['stages'].append(record)
            with (out / (name + '.stdout')).open('w') as stdout, (out / (name + '.stderr')).open('w') as stderr:
                process = subprocess.Popen(record['command'], stdout=stdout, stderr=stderr)
                record['pid'] = process.pid
                save()
                _, status, usage = os.wait4(process.pid, 0)
                process.returncode = os.waitstatus_to_exitcode(status)
            record.update(returncode=process.returncode, wall_seconds=time.perf_counter() - t,
                          user_cpu_seconds=usage.ru_utime, system_cpu_seconds=usage.ru_stime,
                          finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
            save()
            need(process.returncode == 0, name + ' failed; inspect ' + str(out))

        run('field_input', [sys.executable, field / 'verify_tensor_input.py'])
        field_result = json.loads((field / 'independent_verification.json').read_text())
        need(field_result['status'] == 'verified', 'field reconstruction did not verify')
        need(field_result['input_sha256'] == sha(field / 'lattice.json'), 'field reconstruction input binding')
        run('split_certificate', [sys.executable, out / verify_source.name, '--input', field / 'lattice.json',
                                  '--certificate', out / 'certificate.json', '--output', out / 'split_verification'])
        split_result = json.loads((out / 'split_verification/verification.json').read_text())
        need(split_result['status'] == 'verified' and split_result['rational_span_dimension'] == 3916,
             'split certificate did not establish full span')
        for path, digest in report['source_hashes'].items():
            need(sha(Path(path)) == digest, 'source changed while verifying: ' + path)
        need(sha(out / 'certificate.json') == sha(cert), 'certificate copy binding')
        need(sha(field / 'lattice.json') == sha(HERE / 'input/lattice.json'), 'lattice copy binding')
        report.update(status='verified', field_input_verified=True, split_certificate_verified=True,
                      source_and_copy_hashes_unchanged=True, rational_span_dimension=3916,
                      all_seed_norms_8=True, seed_count=split_result['seed_count'],
                      proves_perfection_given_external_minimum_8=True)
    except BaseException as error:
        report.update(status='failed', error=type(error).__name__ + ': ' + str(error))
        raise
    finally:
        report.update(finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      wall_seconds=time.perf_counter() - started,
                      child_user_cpu_seconds=sum(s.get('user_cpu_seconds', 0) for s in report['stages']),
                      child_system_cpu_seconds=sum(s.get('system_cpu_seconds', 0) for s in report['stages']))
        save()
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
