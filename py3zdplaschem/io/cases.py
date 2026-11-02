r"""粒子種・物理設定を含む版管理可能なJSONケースと反応CSVの読み込み。

相対パスは記述元ファイルから解決し、作業ディレクトリに依存させない。JSONには標準ライブラリを使用する。
"""

from __future__ import annotations

from copy import deepcopy
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from ..chemistry.network import NetworkReaction, ReactionNetwork, Species
from ..constants import KB
from ..paths import DATA_DIR
from .sample_inputs import verify_sample_inputs


def read_json(path):
    r"""UTF-8のJSONケースまたはメタデータを読み込む。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。

    Returns
    -------
    dict or list
        JSONを読み込んで得られた値。

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     path = Path(directory) / "example.json"
    ...     _ = path.write_text('{"value": 3}', encoding="utf-8")
    ...     value = read_json(path)
    >>> value
    {'value': 3}
    """
    return json.loads(Path(path).read_text(encoding="utf-8"))


def verify_sources(data, directory):
    r"""参照する元ファイルが記録したハッシュと一致することを検証する。

    Parameters
    ----------
    data : dict
        この処理で使用する構造化された入力メタデータ。
    directory : pathlib.Path
        マニフェスト内の相対パスを解決する基準ディレクトリ。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> import hashlib
    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     base = Path(directory)
    ...     _ = (base / "input.txt").write_bytes(b"example")
    ...     digest = hashlib.sha256(b"example").hexdigest()
    ...     verify_sources(
    ...         {"source_sha256": {"input.txt": digest}}, base
    ...     )
    """
    for relative, expected in data.get("source_sha256", {}).items():
        path = (directory / relative).resolve()
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source hash mismatch: {path}")


@dataclass(frozen=True)
class Case:
    r"""再現可能な計算ケース。条件の変更にはwith_updatesを使用する。
    """

    data: dict
    source: Path

    def with_updates(self, **updates):
        r"""ケースを複製し、指定した設定項目または辞書全体を上書きする。

        Parameters
        ----------
        **updates : dict
            上書き設定。気体・電界・ソルバー・電子の辞書は1階層だけ統合する。

        Returns
        -------
        Case
            気体・電界・ソルバー・電子設定を統合した新しいケース。

        Examples
        --------
        >>> from py3zdplaschem import load_case
        >>> case = load_case("nitrogen")
        >>> changed = case.with_updates(field={"td": 500.0})
        >>> changed.data["field"]["td"], case.data["field"]["td"]
        (500.0, 1000.0)
        """
        data = deepcopy(self.data)
        for key, value in updates.items():
            if key in ("gas", "solver", "field", "electrons"):
                data[key] = {**data.get(key, {}), **deepcopy(value)}
            else:
                data[key] = deepcopy(value)
        return Case(data, self.source)

    def build_network(self):
        r"""ケース入力を検証し、共通の反応ネットワークを構築する。

        Returns
        -------
        ReactionNetwork
            化学量論、単位、初期電荷を検証した反応ネットワーク。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> from py3zdplaschem import load_case
        >>> model = load_case("nitrogen").build_network()
        >>> len(model.species), len(model.reactions)
        (21, 96)
        """
        data = self.data
        allowed = {
            "schema_version",
            "name",
            "mechanism",
            "gas",
            "initial_densities_m3",
            "field",
            "electrons",
            "solver",
            "diagnostics_species",
            "notes",
            "allow_net_charge",
        }
        if data.get("schema_version") != 1 or set(data) - allowed:
            raise ValueError(
                "Unsupported case schema version or unknown case fields."
            )
        gas = data["gas"]
        if set(gas) - {"temperature_K", "pressure_Pa", "mole_fractions"}:
            raise ValueError("Unknown gas settings.")
        tg, pressure = gas["temperature_K"], gas["pressure_Pa"]
        if (
            not np.isfinite(tg)
            or tg <= 0
            or not np.isfinite(pressure)
            or pressure <= 0
        ):
            raise ValueError(
                "Gas temperature and pressure must be positive and finite."
            )
        path = self.source
        mechanism = data["mechanism"]
        if not isinstance(mechanism, dict):
            raise ValueError(
                "Put species metadata and the reaction CSV reference in"
                " mechanism."
            )
        if mechanism.get("schema_version") != 1:
            raise ValueError("Unsupported mechanism schema version.")
        if set(mechanism) - {
            "schema_version",
            "species",
            "reactions",
            "third_bodies",
            "variables",
            "source_sha256",
            "provenance",
            "notes",
        }:
            raise ValueError("Unknown mechanism fields.")
        verify_sources(mechanism, path.parent)
        if mechanism.get("provenance"):
            verify_sample_inputs(
                (path.parent / mechanism["provenance"]).resolve().parent
            )
        species = [Species.from_dict(item) for item in mechanism["species"]]
        names = {s.name for s in species}
        fractions = gas.get("mole_fractions", {})
        if fractions:
            if set(fractions) - names or any(
                not np.isfinite(v) or v < 0 for v in fractions.values()
            ):
                raise ValueError(
                    "Invalid gas mole fractions or unknown background species."
                )
            if not np.isclose(
                sum(fractions.values()), 1.0, rtol=0, atol=1e-12
            ):
                raise ValueError("Gas mole fractions must sum to one.")
            if any(s.charge != 0 for s in species if s.name in fractions):
                raise ValueError(
                    "Background gas fractions must describe neutral species."
                )
        initial = {
            name: value * pressure / (KB * tg)
            for name, value in fractions.items()
        }
        initial.update(data.get("initial_densities_m3", {}))
        reaction_path = (path.parent / mechanism["reactions"]).resolve()
        reactions = []
        with reaction_path.open(encoding="utf-8-sig", newline="") as fp:
            reader = csv.DictReader(fp)
            required = {"id", "equation", "kind", "value", "length_unit"}
            if not reader.fieldnames or not required <= set(reader.fieldnames):
                raise ValueError(
                    f"Reaction CSV requires columns {sorted(required)}."
                )
            optional = {
                "parameters",
                "source_file",
                "line",
                "reference",
                "temperature_min_K",
                "temperature_max_K",
            }
            if set(reader.fieldnames) - required - optional:
                raise ValueError("Unknown reaction CSV columns.")
            for line, row in enumerate(reader, 2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(
                        f"Malformed CSV row {line}: {reaction_path}"
                    )
                row = {key: value.strip() for key, value in row.items()}
                row["line"] = int(row.get("line") or line)
                source = row.get("source_file")
                row["source_file"] = (
                    str((reaction_path.parent / source).resolve())
                    if source
                    else str(reaction_path)
                )
                for key in ("temperature_min_K", "temperature_max_K"):
                    row[key] = float(row[key]) if row.get(key) else None
                reactions.append(NetworkReaction(**row))
        model = ReactionNetwork(
            species,
            reactions,
            temperature_K=tg,
            pressure_Pa=pressure,
            field_td=data.get("field", {}).get("td", 0.0),
            third_bodies=mechanism.get("third_bodies"),
            variables=mechanism.get("variables", ()),
            initial_densities=initial,
            initial_fractions=fractions,
        )
        initial_y = model.initial_state()
        # Explicit override supports intentionally non-neutral test
        # cases.
        charge = abs(model.charge @ initial_y)
        charged = np.abs(model.charge) @ initial_y
        if not data.get("allow_net_charge", False) and charge > max(
            1e-30, 1e-10 * charged
        ):
            raise ValueError(
                "Initial charge is not balanced; set charged densities"
                " explicitly."
            )
        unknown = set(data.get("diagnostics_species", [])) - names
        if unknown:
            raise ValueError(f"Unknown diagnostic species: {unknown}")
        return model

    def electron_settings(self):
        r"""選択した電子バックエンドの設定を解決して検証する。

        Returns
        -------
        dict
            断面積の絶対パスと衝突対応情報を含むバックエンド設定。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> from py3zdplaschem import load_case
        >>> settings = load_case("argon").electron_settings()
        >>> settings["profile_data"]["targets"]
        {'Ar': 'AR'}
        """
        settings = deepcopy(self.data.get("electrons", {"backend": "none"}))
        allowed = {
            "profile",
            "backend",
            "update_mode",
            "update_ns",
            "cells",
            "rtol",
            "max_iterations",
            "mean_energy_eV",
            "coefficients_SI",
            "cross_sections",
        }
        if set(settings) - allowed:
            raise ValueError(
                f"Unknown electron settings: {set(settings) - allowed}"
            )
        backend = settings.setdefault(
            "backend", "bolos" if "profile" in settings else "none"
        )
        if backend not in ("bolos", "none", "prescribed"):
            raise ValueError(
                "Electron backend must be bolos, prescribed or none."
            )
        if backend == "bolos":
            profile_path = self.source
            profile = settings["profile"]
            if not isinstance(profile, dict):
                raise ValueError(
                    "Put collision targets and channel mappings in"
                    " electrons.profile."
                )
            if set(profile) - {
                "cross_sections",
                "targets",
                "channels",
                "source_sha256",
                "notes",
                "coarse_max_eV",
                "adaptive_factor",
            }:
                raise ValueError("Unknown electron profile fields.")
            verify_sources(profile, profile_path.parent)
            paths = profile["cross_sections"]
            profile["cross_sections"] = [
                str((profile_path.parent / p).resolve())
                for p in ([paths] if isinstance(paths, str) else paths)
            ]
            if "cross_sections" in settings:
                paths = settings["cross_sections"]
                profile["cross_sections"] = [
                    str((self.source.parent / p).resolve())
                    for p in ([paths] if isinstance(paths, str) else paths)
                ]
            settings["profile_data"] = profile
            settings["profile_path"] = str(profile_path)
        return settings

    def resolved_data(self):
        r"""入力パスを絶対参照にした、再実行用のケースデータを返す。

        Returns
        -------
        dict
            入力参照先を絶対パスにした独立したケースデータ。

        Examples
        --------
        >>> from pathlib import Path
        >>> from py3zdplaschem import load_case
        >>> data = load_case("argon").resolved_data()
        >>> Path(data["mechanism"]["reactions"]).is_absolute()
        True
        """
        data = deepcopy(self.data)
        base = self.source.parent
        mechanism = data["mechanism"]
        mechanism["reactions"] = str((base / mechanism["reactions"]).resolve())
        if mechanism.get("provenance"):
            mechanism["provenance"] = str(
                (base / mechanism["provenance"]).resolve()
            )
        profile = data.get("electrons", {}).get("profile")
        for block in (mechanism, profile):
            if block and block.get("source_sha256"):
                block["source_sha256"] = {
                    str((base / p).resolve()): value
                    for p, value in block["source_sha256"].items()
                }
        if profile:
            paths = profile["cross_sections"]
            profile["cross_sections"] = [
                str((base / p).resolve())
                for p in ([paths] if isinstance(paths, str) else paths)
            ]
        if data.get("electrons", {}).get("cross_sections"):
            paths = data["electrons"]["cross_sections"]
            data["electrons"]["cross_sections"] = [
                str((base / p).resolve())
                for p in ([paths] if isinstance(paths, str) else paths)
            ]
        return data


def load_case(path):
    r"""ケースJSONまたは標準ケース名air・nitrogen・argonを読み込む。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。

    Returns
    -------
    Case
        相対パス解決のために元ファイルの場所を保持するケース設定。

    Examples
    --------
    >>> case = load_case("nitrogen")
    >>> case.data["name"]
    'nitrogen'
    """
    path = Path(path)
    if str(path) in ("air", "nitrogen", "argon"):
        path = DATA_DIR / "cases" / f"{path}.json"
    path = path.resolve()
    return Case(read_json(path), path)
