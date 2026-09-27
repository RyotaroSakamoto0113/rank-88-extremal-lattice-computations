#!/usr/bin/env python3
"""Check exact identities of the distributed lattice inputs and certificates."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    binding = json.loads((ROOT / 'verification/input_binding.json').read_text())
    for label in ('L', 'L0'):
        assert digest(ROOT / 'data' / label / 'lattice.json') == binding['inputs'][label]
        assert digest(ROOT / 'certificates' / (label + '_perfection.json')) == binding['certificates'][label]
    claims = json.loads((ROOT / 'theorem_code_map.json').read_text())['claims']
    identifiers = {claim['id'] for claim in claims}
    assert len(identifiers) == len(claims)
    for claim in claims:
        assert claim['status'] == 'VERIFIED'
        for dependency in claim['computational_dependencies']:
            assert dependency in identifiers
    for name in ('full_L', 'J0', 'groups', 'factor_claims'):
        report = ROOT / binding['reports'][name]
        assert report.is_file(), report
    print(json.dumps({'status': 'verified', 'lattices': ['L', 'L0'],
                      'computational_claims': len(claims),
                      'scope': 'Input and certificate hashes.'}, indent=2))


if __name__ == '__main__':
    main()
