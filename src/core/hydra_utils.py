"""Helper utilities for working with Hydra metadata.

This module centralizes logic that derives sweep/run specific identifiers
so that file names created during multi-run executions remain informative
and collision free.
"""

from __future__ import annotations

import re
from typing import Iterable, List, Optional

from omegaconf import DictConfig, OmegaConf
from omegaconf.errors import MissingMandatoryValue

try:
    from hydra.core.hydra_config import HydraConfig  # type: ignore
except Exception:  # pragma: no cover - hydra unavailable in some contexts
    HydraConfig = None  # type: ignore

_SANITIZE_PATTERN = re.compile(r"[^A-Za-z0-9._-]+")


def build_run_suffix(
    cfg: DictConfig,
    extra_parts: Optional[Iterable[Optional[str]]] = None,
    include_job_identifier: bool = False,
) -> str:
    """Construct a filesystem-safe suffix describing the current Hydra run.

    Parameters
    ----------
    cfg:
        The Hydra configuration passed to the CLI entry point.
    extra_parts:
        Optional iterable of additional strings (for example a resolved seed
        value) that should be included in the suffix. ``None`` or empty strings
        are ignored.
    include_job_identifier:
        When ``True`` append the Hydra job identifier (``job.id`` or
        ``job.num``) which can be helpful when duplicate override combinations
        appear in a sweep.

    Returns
    -------
    str
        A sanitized identifier that can be appended to file names. If no
        information is available, the string ``"run"`` is returned.
    """

    parts: List[str] = []

    def _append(value: Optional[str]) -> None:
        if value is None:
            return
        text = str(value).strip()
        if not text:
            return
        if text not in parts:
            parts.append(text)

    override_dir = OmegaConf.select(cfg, "hydra.job.override_dirname", default=None)

    hydra_cfg = None
    if HydraConfig is not None:
        try:
            hydra_cfg = HydraConfig.get()
        except (ValueError, MissingMandatoryValue):
            hydra_cfg = None

    if override_dir is None and hydra_cfg is not None:
        override_dir = getattr(hydra_cfg.job, "override_dirname", None)
        if not override_dir:
            overrides = getattr(getattr(hydra_cfg, "overrides", None), "task", None)
            if overrides:
                override_dir = ",".join(overrides)

    _append(override_dir)

    if extra_parts:
        for part in extra_parts:
            _append(part)

    job_id = None
    job_num = None

    if include_job_identifier:
        job_id = OmegaConf.select(cfg, "hydra.job.id", default=None)
        job_num = OmegaConf.select(cfg, "hydra.job.num", default=None)

    if hydra_cfg is not None:
        if job_id is None:
            job_id = getattr(hydra_cfg.job, "id", None)
        if job_num is None:
            job_num = getattr(hydra_cfg.job, "num", None)

    if include_job_identifier:
        if job_id is not None:
            _append(f"job{job_id}")
        elif job_num is not None:
            _append(f"job{job_num}")

    if not parts:
        if hydra_cfg is not None:
            sweep_dir = getattr(getattr(hydra_cfg, "sweep", None), "dir", None)
            if sweep_dir:
                fallback_id = getattr(hydra_cfg.job, "id", None)
                if fallback_id is None:
                    fallback_id = getattr(hydra_cfg.job, "num", None)
                if fallback_id is not None:
                    parts.append(f"job{fallback_id}")

    if not parts:
        return "run"

    raw_suffix = "_".join(parts)
    sanitized = _SANITIZE_PATTERN.sub("_", raw_suffix).strip("_")

    return sanitized or "run"
