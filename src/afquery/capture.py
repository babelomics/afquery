import bisect
import itertools
import pickle
import warnings
import pyranges as pr
import pandas as pd
from .constants import ALL_CHROMS, normalize_chrom
from .models import Technology


class CaptureIndex:
    """Interval index over a capture BED, queried by 1-based position.

    Chromosome keys in ``_index`` are always normalized (see
    :func:`~afquery.constants.normalize_chrom`), because queries reach
    :meth:`covers` already normalized. BED files in the wild use either
    convention ('1' for GRCh37/hs37d5, 'chr1' for hg38), so normalizing on
    construction is what keeps the two sides in agreement.
    """

    _always_covered: bool
    _pr: "pr.PyRanges | None"
    # Per-chrom sorted interval index, keyed by normalized chrom:
    #   {chrom: (sorted_starts, corresponding_ends, running_max_of_ends)}
    # The running max lets covers() answer in O(log n) — see there.
    _index: "dict[str, tuple[list[int], list[int], list[int]]]"

    def __init__(self, *, always_covered: bool = False, pyranges_obj=None):
        self._always_covered = always_covered
        self._pr = pyranges_obj
        self._index = {}
        if pyranges_obj is not None:
            df = pyranges_obj.df
            by_chrom: dict[str, list[tuple[int, int]]] = {}
            for chrom, group in df.groupby("Chromosome", observed=True):
                # Normalizing can merge two groups ('1' and 'chr1') into one key,
                # so accumulate before sorting rather than assigning per group.
                key = normalize_chrom(str(chrom))
                by_chrom.setdefault(key, []).extend(
                    zip(group["Start"].tolist(), group["End"].tolist())
                )
            self._index = {
                chrom: self._pack(pairs) for chrom, pairs in by_chrom.items()
            }

    @staticmethod
    def _pack(pairs: list[tuple[int, int]]) -> tuple[list[int], list[int], list[int]]:
        """Sort (start, end) pairs by start and precompute the running max of ends."""
        pairs = sorted(pairs)
        starts = [s for s, _ in pairs]
        ends = [e for _, e in pairs]
        return starts, ends, list(itertools.accumulate(ends, max))

    @classmethod
    def wgs(cls) -> "CaptureIndex":
        return cls(always_covered=True)

    @classmethod
    def from_bed(cls, bed_path: str) -> "CaptureIndex":
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
            try:
                ranges = pr.read_bed(bed_path)
            except (IndexError, AssertionError, pd.errors.EmptyDataError):
                # pyranges raises rather than returning an empty frame when a BED
                # holds no intervals: IndexError for a zero-byte file, AssertionError
                # for a comment-only one. An empty index is the honest representation
                # — the caller reports it (see describe_capture_problem) instead of
                # the build dying on an opaque traceback.
                return cls()
        return cls(pyranges_obj=ranges)

    def covers(self, chrom: str, pos: int) -> bool:
        """Return True if 1-based pos is within any region for this chrom.

        BED is 0-based half-open: [Start, End) → 1-based: Start < pos <= End
        """
        if self._always_covered:
            return True
        entry = self._index.get(chrom)
        if entry is None:
            return False
        starts, _ends, max_ends = entry
        # Find rightmost interval whose Start < pos (0-based start, so Start < pos_1based)
        # bisect_left gives insertion point for pos in starts; all starts[:idx] < pos
        idx = bisect.bisect_left(starts, pos) - 1
        # Intervals may overlap, so any of starts[:idx+1] could reach pos. max_ends[idx]
        # is the largest End among them, so it alone decides — no backwards scan needed.
        return idx >= 0 and max_ends[idx] >= pos

    @property
    def is_always_covered(self) -> bool:
        """True for the WGS sentinel: every position counts as covered."""
        return self._always_covered

    def indexed_chroms(self) -> set[str]:
        """Normalized chrom names present in the index, unknown contigs included."""
        return set(self._index)

    def known_chroms(self) -> set[str]:
        """Indexed chrom names that are real chromosomes (subset of ALL_CHROMS)."""
        return self.indexed_chroms() & set(ALL_CHROMS)

    def is_empty(self) -> bool:
        """True when the index holds no intervals at all.

        The WGS sentinel is not empty — it covers everything without an index.
        """
        return not self._always_covered and not self._index

    def __getstate__(self) -> dict:
        """Pickle without _pr — nothing outside this class reads the PyRanges object.

        Keeping it put a pandas DataFrame in every WES pickle and made loading depend
        on the pyranges/pandas versions that wrote it. Worse, __init__ normalizes the
        chrom keys of _index but not of _pr, so a pickled _pr is a stale copy of the
        raw BED names — exactly what a future reader would reintroduce this bug from.
        _index alone answers covers() and is enough to rebuild after a normalization
        change.
        """
        state = self.__dict__.copy()
        state.pop("_pr", None)
        return state

    def __setstate__(self, state: dict) -> None:
        """Migrate pickles written by older versions.

        Three generations are repaired here, so existing databases self-heal on
        load with no rebuild:
          1. no _index at all (pre binary-search refactor) — rebuilt from _pr;
          2. un-normalized chrom keys ('1' instead of 'chr1'), which made covers()
             miss every WES technology and silently drop them from AN;
          3. 2-tuple entries lacking the running max of ends.
        """
        self.__dict__.update(state)
        # Pickles written since __getstate__ carry no _pr; keep the attribute present
        # so it stays safe to read on any restored object.
        self.__dict__.setdefault("_pr", None)
        index = getattr(self, "_index", None)

        if index is None:
            self._index = {}
            if getattr(self, "_pr", None) is not None:
                CaptureIndex.__init__(
                    self,
                    always_covered=self._always_covered,
                    pyranges_obj=self._pr,
                )
            return

        needs_migration = any(
            chrom != normalize_chrom(chrom) or len(entry) < 3
            for chrom, entry in index.items()
        )
        if not needs_migration:
            return

        # Rebuild from the stored intervals rather than from _pr: _pr is absent in
        # some pickles, and the intervals themselves are all we need.
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for chrom, entry in index.items():
            starts, ends = entry[0], entry[1]
            by_chrom.setdefault(normalize_chrom(chrom), []).extend(zip(starts, ends))
        self._index = {chrom: self._pack(pairs) for chrom, pairs in by_chrom.items()}

    def save(self, path: str) -> None:
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path: str) -> "CaptureIndex":
        with open(path, "rb") as f:
            return pickle.load(f)


def describe_capture_problem(index: CaptureIndex, tech_name: str) -> "str | None":
    """Return a ready-to-emit message if this index cannot cover anything, else None.

    Build time and load time both need this verdict; sharing one implementation is
    what stops them drifting on what counts as broken. The two failures are reported
    separately because they need different fixes:
      - empty: the BED held no intervals (missing, empty or malformed file);
      - unknown contigs: intervals exist but none is on a real chromosome.
    """
    if index.is_always_covered or index.known_chroms():
        return None
    consequence = (
        "its samples will be counted as uncovered at every position, "
        "lowering AN and inflating AF."
    )
    if index.is_empty():
        return f"Capture regions for technology {tech_name!r} are empty — {consequence}"
    found = sorted(index.indexed_chroms())
    shown = ", ".join(found[:5]) + (", ..." if len(found) > 5 else "")
    return (
        f"Capture regions for technology {tech_name!r} match no known chromosome "
        f"(found: {shown}) — {consequence}"
    )


def load_capture_indices(
    technologies: list[Technology], capture_dir: str
) -> dict[int, "CaptureIndex"]:
    result = {}
    for tech in technologies:
        path = f"{capture_dir}/tech_{tech.tech_id}.pickle"
        try:
            result[tech.tech_id] = CaptureIndex.load(path)
        except FileNotFoundError as exc:
            # Defaulting to WGS here would silently inflate AN, so this must stay fatal
            # — only the message improves.
            raise FileNotFoundError(
                f"Missing capture index for technology {tech.tech_name!r} at {path} — "
                "the database is incomplete; rebuild it or re-run 'afquery update-db'."
            ) from exc
    return result
