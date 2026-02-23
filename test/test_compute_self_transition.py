""" ComputeSelfTransition tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id$"


from openalea.sequence_analysis.data_transform import ComputeSelfTransition
from openalea.sequence_analysis import Sequences
try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data

import pytest

@pytest.fixture
def data():
    return Sequences(get_shared_data("data/sequences1.seq"))

@pytest.fixture
def myi(data):
    # init expect the 4th argument to be provided.
    # vectors is therefore passed as dummy structure
    return interface(data, get_shared_data("data/sequences1.seq"), Sequences)

@pytest.fixture
def datan():
    return Sequences(get_shared_data("data/sequences2.seq"))

@pytest.fixture
def myin(datan):
    # init expect the 4th argument to be provided.
    # vectors is therefore passed as dummy structure
    return interface(data, get_shared_data("data/sequences2.seq"), Sequences)


def test_ComputeSelfTransition():
    ComputeSelfTransition(myi())
    ComputeSelfTransition(myin())

def test_ComputeSelfTransition_order():
    """not implemented see export_markovian_sequences code
    the order arguments is protected..."""
    ComputeSelfTransition(datan(), Order=2)

if __name__ == "__main__":
    test_ComputeSelfTransition()
    test_ComputeSelfTransition_order()
