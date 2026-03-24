import doctest

import matplotlib
import pytest

matplotlib.use("Agg")

import cabaret.camera
import cabaret.focuser
import cabaret.observatory
import cabaret.queries
import cabaret.site
import cabaret.sources
import cabaret.telescope

from .utils import has_internet

pytestmark = pytest.mark.filterwarnings("ignore:.*FigureCanvasAgg.*:UserWarning")


class _LenientOutputChecker(doctest.OutputChecker):
    """Accept unexpected output when none was expected (e.g. Gaia server warnings)."""

    def check_output(self, want, got, optionflags):
        if not want.strip():
            return True
        return super().check_output(want, got, optionflags)


@pytest.mark.parametrize(
    "mod",
    [
        pytest.param(cabaret.sources, id=cabaret.sources.__name__),
        pytest.param(
            cabaret.queries,
            marks=pytest.mark.skipif(not has_internet(), reason="Requires internet"),
            id=cabaret.queries.__name__,
        ),
        pytest.param(cabaret.camera, id=cabaret.camera.__name__),
        pytest.param(cabaret.site, id=cabaret.site.__name__),
        pytest.param(cabaret.focuser, id=cabaret.focuser.__name__),
        pytest.param(cabaret.telescope, id=cabaret.telescope.__name__),
        pytest.param(
            cabaret.observatory,
            marks=pytest.mark.skipif(not has_internet(), reason="Requires internet"),
            id=cabaret.observatory.__name__,
        ),
    ],
)
def test_doctests(mod):
    flags = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
    checker = _LenientOutputChecker()
    runner = doctest.DocTestRunner(optionflags=flags, checker=checker)
    for test in doctest.DocTestFinder().find(mod):
        runner.run(test)
    summary = runner.summarize()
    assert summary.failed == 0, f"{summary.failed} doctest failures in {mod.__name__}"
