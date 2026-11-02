r"""NumPy形式の関数docstringと実行可能な使用例を検査する。
"""

import argparse
import ast
import doctest
from importlib import import_module
import json
import os
from pathlib import Path
import pkgutil
import re
import sys


MODEL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODEL))
source_files = import_module("verification.check_style").source_files


def returns_value(node):
    r"""内部関数の本体を除外し、値を返すreturn文を検出する。

    Parameters
    ----------
    node : ast.AST
        再帰的に調べる関数本体または構文木ノード。

    Returns
    -------
    bool
        関数本体がNone以外の値を返す場合はTrue。
    """
    for child in ast.iter_child_nodes(node):
        if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if isinstance(child, ast.Return):
            if child.value is not None and not (
                isinstance(child.value, ast.Constant)
                and child.value.value is None
            ):
                return True
        if returns_value(child):
            return True
    return False


def check_source(source, label):
    r"""Pythonコード内の各関数のdocstringを検査する。

    Parameters
    ----------
    source : str
        内部関数を含むPythonソースコード。
    label : str
        診断に表示するファイル名と、必要に応じたセル番号。

    Returns
    -------
    tuple
        関数数、Examplesのある関数数、検査メッセージ。
    """
    count, examples, errors = 0, 0, []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        count += 1
        doc = ast.get_docstring(node) or ""
        location = f"{label}:{node.lineno}:{node.name}"
        if not doc:
            errors.append(f"{location}: missing docstring")
            continue
        if not re.search(r"[\u3040-\u30ff\u3400-\u9fff]", doc.splitlines()[0]):
            errors.append(f"{location}: summary must be in Japanese")
        arguments = [
            item.arg
            for item in (
                node.args.posonlyargs
                + node.args.args
                + node.args.kwonlyargs
            )
            if item.arg not in ("self", "cls")
        ]
        if node.args.vararg:
            arguments.append(node.args.vararg.arg)
        if node.args.kwarg:
            arguments.append(node.args.kwarg.arg)
        if arguments:
            section = re.search(
                r"(?ms)^Parameters\n-+\n(.*?)(?=^\w[^\n]*\n-+\n|\Z)",
                doc,
            )
            if section is None:
                errors.append(f"{location}: missing NumPy Parameters")
            else:
                declared = set()
                for header in re.findall(
                    r"(?m)^([\w*, ]+)\s+:\s+\S.*$", section[1]
                ):
                    declared.update(
                        name.strip().lstrip("*")
                        for name in header.split(",")
                    )
                missing = set(arguments) - declared
                extra = declared - set(arguments)
                if missing or extra:
                    errors.append(
                        f"{location}: parameter mismatch: "
                        f"missing={sorted(missing)}, extra={sorted(extra)}"
                    )
        if returns_value(node) and not re.search(
            r"(?m)^(Returns|Yields)\n-+\n", doc
        ):
            errors.append(f"{location}: missing NumPy Returns")
        if re.search(r"(?m)^Examples\n-+\n", doc):
            examples += 1
        try:
            doctest.DocTestParser().get_examples(doc)
        except ValueError as error:
            errors.append(f"{location}: malformed example: {error}")
    return count, examples, errors


def main():
    r"""プロジェクトの関数docstringと指定された使用例を検証する。

    Examples
    --------
    プロジェクトのルートから実行します::

        python verification/check_docstrings.py --examples
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", action="store_true")
    args = parser.parse_args()
    count, examples, errors = 0, 0, []
    for path in source_files(".py"):
        result = check_source(path.read_text(encoding="utf-8-sig"), str(path))
        count += result[0]
        examples += result[1]
        errors.extend(result[2])
    for path in source_files(".ipynb"):
        notebook = json.loads(path.read_text(encoding="utf-8"))
        for index, cell in enumerate(notebook["cells"], 1):
            if cell["cell_type"] == "code":
                result = check_source(
                    "".join(cell["source"]), f"{path}:cell-{index}"
                )
                count += result[0]
                examples += result[1]
                errors.extend(result[2])
    attempted, failed = 0, 0
    if args.examples and not errors:
        os.environ.setdefault("MPLBACKEND", "Agg")
        package = import_module("plasglo")
        for item in pkgutil.walk_packages(
            package.__path__, package.__name__ + "."
        ):
            if item.name.endswith(".__main__"):
                continue
            result = doctest.testmod(
                import_module(item.name), optionflags=doctest.ELLIPSIS
            )
            attempted += result.attempted
            failed += result.failed
    report = dict(
        functions=count,
        functions_with_examples=examples,
        documentation_errors=errors,
        executable_examples_attempted=attempted,
        executable_examples_failed=failed,
    )
    output = MODEL / "outputs/verification/pep8/docstrings.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    for error in errors:
        print(error)
    print(
        f"{count} function docstrings, {examples} with Examples; "
        f"{len(errors)} documentation errors; "
        f"{failed}/{attempted} executable examples failed."
    )
    raise SystemExit(bool(errors or failed))


if __name__ == "__main__":
    main()
