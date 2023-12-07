from __future__ import annotations

import importlib.metadata

import asymmetric_immiscibility_simulations as m


def test_version():
    assert importlib.metadata.version("asymmetric_immiscibility_simulations") == m.__version__
