import json
from dataclasses import asdict, is_dataclass
from typing import Any, Union, TextIO
import pint

# ----- Internals -----
class _PintEncoder(json.JSONEncoder):
    def default(self, obj: Any):
        # pint Quantity -> tagged dict
        if isinstance(obj, pint.Quantity):
            mag = obj.magnitude
            # make numpy arrays/lists JSON-friendly if present
            if hasattr(mag, "tolist"):
                mag = mag.tolist()
            return {
                "__pint_quantity__": True,
                "magnitude": mag,
                "units": str(obj.units),
            }
        # Expand dataclasses to dicts
        if is_dataclass(obj):
            return asdict(obj)
        return super().default(obj)

def _pint_object_hook_factory(ureg: pint.UnitRegistry):
    def hook(d: dict):
        if d.get("__pint_quantity__") is True:
            return ureg.Quantity(d["magnitude"], ureg.parse_units(d["units"]))
        return d
    return hook

# ----- Public API -----
def json_pint_dump(
    obj: Any,
    fp: Union[str, TextIO],
    ureg: pint.UnitRegistry,
    **json_kwargs,
) -> None:
    """
    Save obj to JSON, converting pint.Quantity to a tagged dict.
    - fp: path or open file-like in text mode
    - ureg: your UnitRegistry (not used at dump-time but kept symmetric)
    Extra json kwargs are passed to json.dump (indent, sort_keys, etc.).
    """
    # ensure_ascii False by default (override if you pass ensure_ascii=...)
    json_kwargs.setdefault("ensure_ascii", False)
    if isinstance(fp, (str, bytes, bytearray)):
        encoding = json_kwargs.pop("encoding", "utf-8")
        with open(fp, "w", encoding=encoding) as f:
            json.dump(obj, f, cls=_PintEncoder, **json_kwargs)
    else:
        json.dump(obj, fp, cls=_PintEncoder, **json_kwargs)

def json_pint_load(
    fp: Union[str, TextIO],
    ureg: pint.UnitRegistry,
):
    """
    Load JSON and reconstruct pint.Quantity using the provided UnitRegistry.
    - fp: path or open file-like in text mode
    - ureg: your UnitRegistry
    """
    hook = _pint_object_hook_factory(ureg)
    if isinstance(fp, (str, bytes, bytearray)):
        with open(fp, "r", encoding="utf-8") as f:
            return json.load(f, object_hook=hook)
    else:
        return json.load(fp, object_hook=hook)
