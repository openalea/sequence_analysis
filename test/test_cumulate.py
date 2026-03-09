"""Cumulate tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__revision__ = "$Id$"


from openalea.sequence_analysis.data_transform import Cumulate
from .tools import robust_path as get_shared_data
from openalea.sequence_analysis import Sequences

seqn = Sequences(str(get_shared_data("sequences2.seq")))
seq1 = Sequences(str(get_shared_data("sequences1.seq")))


def test_cumulate1():
    data = seq1
    a = Cumulate(data)
    b = data.cumulate(1).markovian_sequences()

    assert str(a) == str(b)
    assert a.get_max_value(0) == 29
    assert b.get_max_value(0) == 29

    assert data.get_max_value(0) == 2


def test_cumulaten():
    for var in range(1, seqn.nb_variable + 1):
        assert str(seqn.cumulate(var).markovian_sequences()) == str(Cumulate(seqn, var))
