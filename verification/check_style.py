r"""プロジェクトのPythonコードとノートブックのPEP 8準拠を検査する。
"""

import json
from pathlib import Path

import pycodestyle

MODEL = Path(__file__).resolve().parents[1]
EXCLUDED = {"Py3bolos", "outputs", ".git", "__pycache__"}


def source_files(suffix):
    r"""保存コードとルートのスクリプトを含む、管理対象のソースを選ぶ。

    Parameters
    ----------
    suffix : str
        先頭のドットを含む、対象ソースファイルの拡張子。

    Returns
    -------
    list of pathlib.Path
        依存ライブラリと出力を除外した、並べ替え済みソース一覧。

    Examples
    --------
    >>> all(path.suffix == ".py" for path in source_files(".py"))
    True
    """
    paths = [
        path
        for path in MODEL.rglob("*" + suffix)
        if not EXCLUDED.intersection(path.relative_to(MODEL).parts)
    ]
    if suffix == ".py":
        paths.extend(MODEL.parent.glob("*.py"))
    return sorted(paths)


def main():
    # W503 conflicts with PEP 8's preferred breaks before binary
    # operators.
    r"""プロジェクトのPythonファイルとノートブックをPEP 8で検査する。

    Raises
    ------
    SystemExit
        コマンドの終了時、または検査エラーを報告する場合。

    Examples
    --------
    プロジェクトのルートから実行します::

        python verification/check_style.py
    """
    style = pycodestyle.StyleGuide(
        ignore=["W503"],
        max_line_length=79,
        max_doc_length=72,
    )
    sources = source_files(".py")
    for path in sources:
        style.input_file(str(path))
    notebooks = source_files(".ipynb")
    code_cells = 0
    for path in notebooks:
        notebook = json.loads(path.read_text(encoding="utf-8"))
        for index, cell in enumerate(notebook["cells"], 1):
            if cell["cell_type"] != "code":
                continue
            source = "".join(cell["source"]).rstrip() + "\n"
            style.input_file(
                f"{path}:cell-{index}",
                lines=source.splitlines(True),
            )
            code_cells += 1
    count = style.options.report.total_errors
    print(
        f"{len(sources)} Python files, {len(notebooks)} notebooks, "
        f"{code_cells} code cells: {count} style violations."
    )
    raise SystemExit(bool(count))


if __name__ == "__main__":
    main()
