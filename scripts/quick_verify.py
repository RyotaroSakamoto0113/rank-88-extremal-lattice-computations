#!/usr/bin/env python3
"""Check the distribution and explicit certificates; no exhaustive search."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]

def check_hashes():
    manifest = ROOT / 'SHA256SUMS'
    if not manifest.is_file():
        raise RuntimeError('SHA256SUMS is missing')
    count = 0
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        path = ROOT / name
        if not path.resolve().is_relative_to(ROOT):
            raise RuntimeError('Unsafe manifest path')
        h = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(block)
        if h.hexdigest() != digest:
            raise RuntimeError('Hash mismatch: ' + name)
        count += 1
    return count

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, help='New output directory')
    args = ap.parse_args()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='rank88-quick-'))
    if args.output:
        out.mkdir(parents=True, exist_ok=False)
    result = {'scope': 'Hashes, exact formulas, Gram matrices and separate perfection certificates.',
              'files_hashed': check_hashes()}
    binding = subprocess.run([sys.executable, ROOT/'scripts/check_input_binding.py'],
                             capture_output=True, text=True, check=True)
    result['input_binding'] = json.loads(binding.stdout)
    formulas = subprocess.run([sys.executable, ROOT/'scripts/verify_formulas.py'], capture_output=True, text=True, check=True)
    result['formulas'] = json.loads(formulas.stdout)
    for label in ['L', 'L0']:
        subprocess.run([sys.executable, ROOT/'scripts/verify_split_certificate.py',
                        '--input', ROOT/'data'/label/'lattice.json',
                        '--certificate', ROOT/'certificates'/(label+'_perfection.json'),
                        '--output', out/label], check=True, stdout=subprocess.DEVNULL)
        record = json.loads((out/label/'verification.json').read_text())
        assert record['status'] == 'verified'
        result[label] = {k: record[k] for k in ['status','seed_count','rational_span_dimension','positive_even_unimodular_gram_verified']}
    result['status'] = 'quick_verified'
    (out/'quick_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    print('Output:', out)

if __name__ == '__main__':
    main()
