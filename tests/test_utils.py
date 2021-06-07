# python modules
import pytest
# test modules
from gsc.utils import compile_maps


@pytest.mark.parametrize('maps, expected_map', [
    ([{}], {}),
    ([{}, {}], {}),
    ([{0: 1}, {1: 2}], {0: 2}),
    ([{0: 1, 1: 2, 2: 0}, {0: 1, 1: 2, 2: 0}], {0: 2, 1: 0, 2: 1}),
])
def test_compile_maps(maps, expected_map):
    assert compile_maps(maps) == expected_map
