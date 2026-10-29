"""Run report saved next to each output as a text file: what was done, with what, and what was skipped."""
import os
import platform
import sys
from datetime import datetime

from . import __version__


def report_path(output):
    return output + ".report.txt"


def _lines(section, indent=2):
    # align values of plain items; headings of nested sections do not count
    width = max((len(str(key)) for key, value in section.items() if not isinstance(value, dict)), default=0) + 1
    for key, value in section.items():
        if isinstance(value, dict):
            yield f"{' ' * indent}{key}:"
            yield from _lines(value, indent + 2) if value else [f"{' ' * (indent + 2)}none"]
        elif isinstance(value, (list, tuple)):
            yield f"{' ' * indent}{(str(key) + ':').ljust(width)} {', '.join(map(str, value)) or 'none'}"
        else:
            yield f"{' ' * indent}{(str(key) + ':').ljust(width)} {value}"


def write_report(output, title, sections):
    """
    Write <output>.report.txt. sections: {heading: {item: value}}; a value can be
    a nested dict (indented sub-items) or a list.
    """
    header = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
        "NANOSLAST": __version__,
        "Python": platform.python_version(),
        "output": os.path.abspath(output),
    }
    lines = [f"NANOSLAST {title} report", "=" * (len(title) + 17), *_lines(header, 0)]
    for heading, section in sections.items():
        lines += ["", heading, "-" * len(heading), *_lines(section)]

    with open(report_path(output), "w") as f:
        f.write("\n".join(lines) + "\n")
    return report_path(output)
