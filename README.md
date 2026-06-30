# Sanger Primer Designer

環状プラスミドのサンガーシークエンス確認用プライマーセットを設計するローカルGUIツールです。

FASTAまたはGenBank形式のプラスミド配列を入力し、指定したdepthを満たすように既存プライマーを選択し、不足する場合は新規プライマーを設計します。

## 必要なもの

- Python 3.11以上
- `uv`

依存関係は `uv` で管理します。

```powershell
uv sync
```

## GUIの起動

リポジトリのルートディレクトリで以下を実行します。

```powershell
uv run streamlit run src/sanger_designer/gui.py
```

ブラウザで以下を開きます。

```text
http://localhost:8501
```

終了する場合は、Streamlitを起動したターミナルで `Ctrl+C` を押します。

## 入力ファイル

### プラスミド配列

GUIの `Plasmid sequence` に、以下のいずれかをアップロードします。

- GenBank: `.gb`, `.gbk`
- FASTA: `.fa`, `.fasta`
- 配列テキスト: `.txt`

入力配列は常に環状プラスミドとして扱います。GenBankのfeature注釈は使わず、配列のみを読み込みます。

FASTAとGenBankは単一レコードのみ対応します。FASTAで複数のヘッダ行がある場合、またはGenBankで `ORIGIN` が複数ある場合はエラーになります。GenBankで `ORIGIN` 内に配列がない場合もエラーになります。

### 既存プライマーリスト

GUIの `Existing primer lists` に、既存プライマーリストをアップロードできます。複数ファイルを同時に指定できます。

`Use default_primers.txt` を有効にすると、ローカルの標準プライマーリストも併用します。

```text
src/sanger_designer/default_primers.txt
```

このファイルはユーザーが編集するローカル設定ファイルです。Git管理には含めません。存在しない場合、GUIは警告を表示し、既存プライマーなしで処理を続行します。

標準リストとアップロードリストを併用した場合、同じ配列のプライマーが複数存在すると標準リスト側の名前を優先します。アップロードリスト同士では、先に読み込まれたリストの名前を優先します。

## プライマーリスト形式

以下のいずれかの区切り形式を使えます。

```text
primer_name<TAB>sequence<TAB>memo
primer_name,sequence,memo
primer_name;sequence;memo
```

`memo` は省略できます。

例:

```text
MyPrimer_F	ACGTACGTACGTACGTACGT	optional memo
MyPrimer_R	TGCATGCATGCATGCATGCA
```

## 基本的な使い方

1. `Plasmid sequence` に配列ファイルをアップロードします。
2. 必要に応じて `Existing primer lists` に既存プライマーリストをアップロードします。
3. 標準リストを使う場合は `Use default_primers.txt` を有効にします。
4. `Target depth` を選択します。通常はデフォルトの `2` を使います。
5. プライマーを置きたくない領域がある場合は `Mask regions` に入力します。
6. 新規プライマー名のprefixを変える場合は `New primer prefix` を変更します。
7. `Design primers` を押します。

## マスク指定

マスク領域にはプライマーを結合させません。ただし、マスク領域もリードで読まれる必要があります。

形式:

```text
63..88
63..88,3025..3063
```

環状配列のoriginを跨ぐ指定もできます。

```text
7200..120
```

## 詳細設定

`Advanced settings` では以下を変更できます。

- Read length
- Noisy bases from primer
- Minimum binding gap
- Preferred pair min/max
- Primer length min/max
- Tm min/max
- GC min/max

通常はデフォルト値のままで使えます。

## 結果の見方

### Overview

設計が成功したか、最低depth、使用プライマー数を確認できます。プライマーリストとdepthレポートを `.txt` としてダウンロードできます。

### Depth

各塩基位置のdepthを線グラフで表示します。赤い領域があれば、指定depthを満たしていない不足領域です。

### Placement

プライマーの結合位置とリードで読まれる範囲を直線上に表示します。

- 右向き三角: Forward
- 左向き三角: Reverse
- 点線: プライマー近傍のノイズ領域
- 太いバー: 採用されるリード範囲

### Primers

採用されたプライマーを表で確認できます。既存プライマーと新規設計プライマーで絞り込みできます。

### Report

結果概要と実行時パラメータをテキストで確認できます。

## サポートバンドル

IssueやPull Requestで再現情報を共有するために、`Overview` からサポートバンドルをzip形式でダウンロードできます。

Bundleにはプラスミド配列やプライマー配列が入ります。共有前に内容を確認してください。

バンドルには以下の情報が含まれます。

- 実行時の設定値
- `defaults.py` のデフォルト値
- `default_primers.txt` の内容
- アップロードしたプラスミド配列
- アップロードした既存プライマーリスト
- マージ後の既存プライマーリスト
- 採用プライマーリスト
- depthレポート
- positionごとのdepthCSV
- 実行ログと警告

## 出力

出力プライマーリストは以下の情報を含みます。

- プライマー名
- 配列
- メモ欄: position、direction、binding、cover

最終出力には、Tm値、GC含有率、既存/新規の区別は含めません。

## トラブルシュート

### `default_primers.txt` が見つからない

`Use default_primers.txt` が有効で、以下のファイルが存在しない場合に警告が出ます。

```text
src/sanger_designer/default_primers.txt
```

標準リストを使いたい場合は、この場所にプライマーリストを作成してください。標準リストを使わない場合は、チェックを外すか、`Existing primer lists` に任意のリストをアップロードしてください。

### depthを達成できない

マスク、ユニーク性、Tm/GC、プライマー間隔の制約により、指定depthを達成できないことがあります。その場合でも、ツールは最大限の候補セットと不足領域を返します。

### 起動できない

依存関係が不足している可能性があります。

```powershell
uv sync
uv run streamlit run src/sanger_designer/gui.py
```

## 現在の制限

- 入力配列は環状プラスミドとして扱います。
- GenBank feature注釈は使いません。
- ユニーク性判定は完全一致ベースです。
- 曖昧塩基は想定していません。

