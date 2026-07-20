import logging

from ..capture import CaptureIndex
from ..constants import ALL_CHROMS
from ..models import Technology

logger = logging.getLogger(__name__)


def build_capture_indices(
    technologies: list[Technology],
    capture_dir: str,
) -> None:
    for tech in technologies:
        if tech.bed_path is None:
            idx = CaptureIndex.wgs()
        else:
            idx = CaptureIndex.from_bed(tech.bed_path)
            if not any(c in ALL_CHROMS for c in idx._index):
                logger.warning(
                    "[regions] BED for technology '%s' (%s) matches no known chromosome — "
                    "its samples would be counted as uncovered at every position.",
                    tech.tech_name, tech.bed_path,
                )
        idx.save(f"{capture_dir}/tech_{tech.tech_id}.pickle")
