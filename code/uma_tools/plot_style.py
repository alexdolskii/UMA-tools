"""Shared typography and export settings for UMA scientific figures."""

from functools import lru_cache

FONT_SIZES = {"title": 20, "panel": 16, "axis": 14, "sample": 12, "note": 12}
PNG_DPI = 300
HEADER_COLOR = "234F68"


@lru_cache(maxsize=1)
def plot_font():
    """Use the FIA font preference and a portable fallback."""
    from matplotlib.font_manager import fontManager

    installed = {font.name for font in fontManager.ttflist}
    return next(
        font
        for font in ("Arial", "Liberation Sans", "DejaVu Sans")
        if font in installed
    )


def rc_parameters():
    """Keep labels literal and embed TrueType fonts in vector PDF."""
    return {
        "font.family": plot_font(),
        "font.size": FONT_SIZES["axis"],
        "pdf.fonttype": 42,
        "text.usetex": False,
        "text.parse_math": False,
    }


@lru_cache(maxsize=2048)
def wrap_label(text, width_points, size, weight="normal"):
    """Wrap to measured width without dropping oversized tokens."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    font = FontProperties(family=plot_font(), size=size, weight=weight)

    def width(value):
        return TextPath((0, 0), value, prop=font).get_extents().width

    lines = []
    for paragraph in text.split("\n"):
        line = ""
        for word in paragraph.split():
            candidate = f"{line} {word}".strip()
            if line and width(candidate) > width_points:
                lines.append(line)
                line = word
            else:
                line = candidate
            while len(line) > 1 and width(line) > width_points:
                low, high = 1, len(line)
                while low + 1 < high:
                    middle = (low + high) // 2
                    if width(line[:middle]) <= width_points:
                        low = middle
                    else:
                        high = middle
                lines.append(line[:low])
                line = line[low:]
        lines.append(line)
    return "\n".join(lines)


def short_labels(groups):
    """Move shared words to the title, retaining group identities."""
    tokens = [group.split() for group in groups]
    shared = 0
    for words in zip(*tokens):
        if len(set(words)) != 1:
            break
        shared += 1
    shared = min(shared, min(map(len, tokens)) - 1)
    labels = {
        group: " ".join(words[shared:]) for group, words in zip(groups, tokens)
    }
    if shared <= 0 or len(set(labels.values())) != len(groups):
        return "", {group: group for group in groups}
    return " ".join(tokens[0][:shared]), labels


def panel_letter(index):
    """Return stable panel letters, including AA after Z."""
    result = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        result = chr(65 + remainder) + result
    return result
