#!/usr/bin/env python3
"""Hash the external archive and every regular member without extracting."""
from pathlib import Path, PurePosixPath
import argparse,hashlib,json,tarfile
ROOT=Path(__file__).resolve().parents[1]
def digest(stream):
    h=hashlib.sha256()
    for b in iter(lambda:stream.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);a=ap.parse_args()
    manifest=json.loads((ROOT/'zenodo-large-file-manifest.json').read_text())
    record=next(x for x in manifest['assets']if x['filename']==a.archive.name)
    assert a.archive.stat().st_size==record['bytes'],'Archive size mismatch'
    with a.archive.open('rb')as f:assert digest(f)==record['sha256'],'Archive hash mismatch'
    entries={x['path']:x for x in map(json.loads,(ROOT/'full_certificate_index.jsonl').read_text().splitlines())}
    seen=set()
    with tarfile.open(a.archive,'r|gz')as tf:
        for m in tf:
            p=PurePosixPath(m.name)
            assert not p.is_absolute() and '..'not in p.parts,'Unsafe archive path'
            if m.isdir():continue
            assert m.isfile() and m.name in entries and m.name not in seen,'Unexpected member'
            seen.add(m.name);e=entries[m.name]
            assert m.size==e['bytes'],'Member size mismatch: '+m.name
            with tf.extractfile(m)as f:assert digest(f)==e['sha256'],'Member hash mismatch: '+m.name
    assert seen==set(entries),'Missing archive members'
    print(json.dumps({'status':'verified','archive':a.archive.name,'members':len(seen)},indent=2))
if __name__=='__main__':main()
