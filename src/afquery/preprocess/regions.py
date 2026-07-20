import warnings

from ..capture import CaptureIndex, describe_capture_problem
from ..models import AfqueryWarning, Technology


def build_capture_indices(
    technologies: list[Technology],
    capture_dir: str,
) -> None:
    for tech in technologies:
        if tech.bed_path is None:
            idx = CaptureIndex.wgs()
        else:
            idx = CaptureIndex.from_bed(tech.bed_path)
            problem = describe_capture_problem(idx, tech.tech_name)
            if problem is not None:
                # warnings, not logger: this is also reached from the Python API via
                # update.py, where a logger message is invisible unless the caller
                # configured a handler. "Your allele frequencies will be wrong" has to
                # surface on both paths.
                warnings.warn(
                    f"{problem} (BED: {tech.bed_path})", AfqueryWarning, stacklevel=2,
                )
        idx.save(f"{capture_dir}/tech_{tech.tech_id}.pickle")
