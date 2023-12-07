from io import StringIO
import ais.utils
import numpy as np
import pytest


def test_chain_parse_single_valence(EK_paramfile):
    seqfile = StringIO('EKEKEKEK')
    result = chain_parse(seqfile, 1, EK_paramfile)
    assert len(result) == 5
    assert result[0] == [0, 1] * 4
    assert result[1] == pytest.approx([0, 1] * 4)


@pytest.fixture
def EK_paramfile():
    return StringIO(
        '#AA     Mass    Charge  Sigma   Lambda PROLINE MODDED\n'
        'GLU     129.10  -1.00   5.920   0.459\n'
        'LYS     128.20  1.00    6.360   0.514\n'
        '#another comment\n'
    )


@pytest.fixture
def small_paramfile():
    return StringIO(
        '#AA     Mass    Charge  Sigma   Lambda PROLINE MODDED\n'
        'ALA     71.08   0.00    5.040   0.730\n'
        'ARG     156.20  1.00    6.560   0.000\n'
        '#another comment\n'
    )

def test_paramparse(small_paramfile):

    result = paramparse(small_paramfile)
    assert len(result) == 6
    assert result[0] == pytest.approx([71.08, 156.20])
    assert result[1] == pytest.approx([0, 1])
    assert result[2] == pytest.approx([5.040, 6.560])
    assert result[3] == pytest.approx([0.730, 0])
    assert result[4] == ['ALA', 'ARG']

    _param_matches_smallfile(result[5])


def test_get_param_dict(small_paramfile):
    result = get_param_dict(small_paramfile)
    _param_matches_smallfile(result)


def _param_matches_smallfile(result):
    assert len(result) == 2
    assert list(result.keys()) == ['ALA', 'ARG']
    assert result['ALA'] == pytest.approx([71.08, 0.00, 5.040, 0.730])
    assert result['ARG'] == pytest.approx([156.20, 1.00, 6.560, 0.000])
