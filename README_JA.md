# 88次元 extremal even unimodular lattices の計算資料

このリポジトリには、2つの格子 `L` と `L0` の構成に用いた計算コード、厳密な入力データ、検証プログラム、小規模な証明書を収録します。大規模な全列挙の証明書は別アーカイブとし、全ファイルのサイズとSHA-256を `full_certificate_index.jsonl` に記録します。

## 短時間の検算

Linux/macOS の Python 3.9 以上で実行します。

```sh
python3 scripts/quick_verify.py
```

配布ファイルのハッシュ、入力の同一性、厳密な数式計算、Gram行列の正定値性・偶性・行列式1、両格子のノルム8の証人と独立したperfection証明書を確認します。

## 全計算の再実行

PARI/GP、C++17コンパイラ、GMP、GAP（SmallGrpとGAPDoc）を用意し、存在しない出力ディレクトリを指定します。

```sh
python3 scripts/reproduce_all.py --output /absolute/path/to/new-rank88-run
```

数GB以上の空き領域を用意してください。個別の実行方法と必要な環境は [再現手順](docs/reproducibility.md) を参照してください。

## 内容

- `data/`：両格子の厳密なGram行列、作用、因子群入力。
- `certificates/`：両格子を別々に検査できるperfection証明書。
- `compute/`：`L`、`J0`、軌道表、類数、単数群、因子群の計算プログラムと入力。
- `scripts/`：入力の再構成、短時間の検算、全計算の再実行。
- `verification/`：検証結果の小規模な記録。
- `theorem_code_map.*` と `computation_dependency_map.*`：計算主張とコード・入力・証拠・期待出力の対応。
- `full_certificate_index.jsonl`：別アーカイブ内の全ファイルのサイズとSHA-256。
- `zenodo-large-file-manifest.json`：アーカイブの識別情報とハッシュ値。

全計算の新規実行には、大容量アーカイブは不要です。計算上の旧ディレクトリ名は、再実行プログラムとの整合のため保持します。

コードのライセンスは `LICENSE`、第三者資料の権利は `NOTICE.md` を参照してください。
