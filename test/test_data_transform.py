"""various function tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__revision__ = "$Id$"

import pytest

from openalea.sequence_analysis import (
    ComputeStateSequences,
    RemoveRun,
    Sequences,
    TransitionCount,
)
from openalea.stat_tool import (
    Display,
    Distribution,
    Merge,
    Mixture,
    Plot,
    SelectStep,
    Simulate,
    Vectors,
)

from openalea.sequence_analysis.data_transform import (
    Cumulate,
    Difference
)

from openalea.stat_tool.data_transform import *
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.cluster import Transcode, Cluster

try:
    from .tools import DISABLE_PLOT, interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import DISABLE_PLOT, interface
    from tools import robust_path as get_shared_data

@pytest.fixture
def build_data():
    """todo: check identifier output. should be a list"""
    # build a list of 2 sequences with a variable that should be identical
    # to sequences1.seq
    data = Sequences([[1,0,0,0,1,1,2,0,2,2,2,1,1,0,1,0,1,1,1,1,0,1,1,1,0,1,2,2,2,1],
            [0, 0, 0, 1, 1, 0, 2, 0, 2, 2, 2, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0]]
    )
    assert data
    assert data.nb_sequence == 2
    assert data.nb_variable == 1
    assert data.cumul_length == 52
    assert data.max_length == 30
    assert [0, 1] == data.get_identifiers()

    return data


@pytest.fixture
def create_data_sequence():
    return Sequences(str(get_shared_data("sequences1.seq")))


@pytest.fixture
def create_data_sequence2():
    return Sequences(str(get_shared_data("sequences2.seq")))

@pytest.fixture
def build_seqn():
    return Sequences([[[1, 1, 1], [12, 12, 12]], [[2, 2, 2], [22, 23, 24]]])


@pytest.fixture
def build_seq1():
    return Sequences([[1, 1, 1], [2, 2, 2]])


@pytest.fixture
def build_seq_realn():
    return Sequences(
        [
            [[1.5, 1.5, 1.5], [12.5, 12.5, 12.5]],
            [[2.5, 2.5, 2.5], [22.5, 23.5, 24.5]],
        ]
    )


class TestRemoveRun:
    """
    def test_sequences_1(self, create_data_sequence):
        seq1 = create_data_sequence
        seq2 = seq1.remove_run(1, 0, "e", 2)
        seq3 = RemoveRun(seq1, 1, 0, "e", MaxLength=2)
        assert str(seq3) == str(seq2)
    """

    def test_incorrect_value(self, create_data_sequence):
        seq1 = create_data_sequence
        try:
            seq1.remove_run(1, 3, "e", 10)
            assert False
        except:
            assert True
        try:
            seq1.remove_run(1, -1, "e", 10)
            assert False
        except:
            assert True

    def test_incorrect_variable(self, create_data_sequence):
        seq1 = create_data_sequence
        try:
            seq1.remove_run(0, 2, "e", 10)
            assert False
        except:
            assert True
    """
    def test_sequences_2(self, create_data_sequence2):
        seq1 = create_data_sequence2
        seq2 = seq1.remove_run(1, 0, "e", 2)
        seq3 = RemoveRun(seq1, 1, 0, "e", 2)
        assert str(seq3) == str(seq2)
    """

def test_markov_data():
    """not implemented"""
    pass


def test_semi_markov_data():
    """not implemented"""
    pass


def test_discrete_sequences():
    """not implemented"""
    pass


def test_compute_state_sequence():
    from openalea.sequence_analysis import HiddenSemiMarkov

    seq = Sequences(str(get_shared_data("wij1.seq")))
    hsmc0 = HiddenSemiMarkov(str(get_shared_data("wij1.hsc")))
    ComputeStateSequences(seq, hsmc0, Algorithm="ForwardBackward", Characteristics=True)

"""
def test_transition_count():
    seq = Sequences(str(get_shared_data("wij1.seq")))
    TransitionCount(seq, 5, Begin=True, Estimator="MaximumLikelihood", Filename="ASCII")
"""

def test_merge():
    mixt1 = Mixture(
        0.6, Distribution("B", 2, 18, 0.5), 0.4, Distribution("NB", 10, 10, 0.5)
    )

    mixt_histo1 = Simulate(mixt1, 200)

    histo10 = mixt_histo1.extract_component(1)
    histo11 = mixt_histo1.extract_component(2)

    histo12 = Merge(histo10, histo11)

    assert histo12


def test_value_select(build_seqn):
    "test_value_select implemented but need to be checked"
    seqn = build_seqn
    a = seqn.value_select(1, 1, 2, True)
    assert a
    assert str(ValueSelect(seqn, 1, 1, 2)) == str(seqn.value_select(1, 1, 2, True))


def test_select_variable_int(build_seqn):
    "test_select_variable_int implemented but need to be checked (index issue)"
    # !!!!!!!! NEED to CHECK THE INDEX 0, 1 , ... or 1,2,....
    # what about identifiers ?
    # Variable seems to start at 1 not 0
    s = build_seqn
    # select variable 1
    select = s.select_variable([1], keep=True)
    assert select[0, 0] == [1]
    assert select[1, 0] == [2]


def test_select_variable_real(build_seq_realn):
    "test_select_variable_real implemented but need to be checked (index issue)"
    # !!!!!!!! NEED to CHECK THE INDEX 0, 1 , ... or 1,2,....
    # what about identifiers ?
    # Variable seems to start at 1 not 0
    s = build_seq_realn
    # select variable 1
    select = s.select_variable([1], keep=True)
    assert select[0, 0] == [1.5]
    assert select[1, 0] == [2.5]


def test_select_individual(build_seqn):
    # select one or several sequences
    s = build_seqn

    # select all
    select = s.select_individual([0, 1], keep=True)
    assert s.display() == select.markovian_sequences().display()

    # select first sequence only
    select = s.select_individual([0], keep=True)
    assert select[0, 0] == [1, 1, 1]
    assert select[0, 1] == [12, 12, 12]
    try:
        select[1, 0]
        assert False
    except:
        assert True

    # select second sequence only
    select = s.select_individual([1], keep=True)
    assert select[0, 0] == [2, 2, 2]
    assert select[0, 1] == [22, 23, 24]
    try:
        select[1, 0]
        assert False
    except:
        assert True


def test_shift_seqn(build_seqn):
    s = build_seqn
    shifted = s.shift(1, 2)
    assert shifted[0, 0] == [3, 1, 1]
    assert shifted[0, 1] == [14, 12, 12]
    assert shifted[1, 0] == [4, 2, 2]
    assert shifted[1, 1] == [24, 23, 24]
    shifted1 = s.shift(1, 1)
    assert shifted1[0, 0] == [2, 1, 1]
    assert shifted1[0, 1] == [13, 12, 12]
    assert shifted1[1, 0] == [3, 2, 2]
    assert shifted1[1, 1] == [23, 23, 24]


def test_shift_seq1(build_seq1):
    s = build_seq1
    shifted = s.shift(1, 2)
    assert shifted[0, 0] == [3, 1, 1]
    assert shifted[0, 1] == [4, 2, 2]


def test_threshold_seq1(build_seqn):
    s = build_seqn
    thresholded = s.thresholding(1, 10, "ABOVE")
    for x in thresholded:
        for v in x:
            assert v[0] <= 10


def test_threshold_seq():
    s = Sequences(
        [
            [[1.01, 1.07], [2.01, 2.07], [1.99, 1.07], [2.41, 2.07]],
            [[1.97, 1.07], [1.98, 2.07], [1.99, 1.07], [2.00, 2.07]],
        ]
    )
    thresholded = s.thresholding(1, 1.99, "ABOVE")
    for x in thresholded:
        for v in x:
            assert v[0] <= 1.99
    thresholded = s.thresholding(1, 1.99, "BELOW")
    for x in thresholded:
        for v in x:
            assert v[0] >= 1.99

    assert s


def test_merge(build_seqn, build_seq1):
    s1 = build_seqn
    s2 = build_seqn
    s3 = build_seq1

    sall = s1.merge([s2])
    assert sall.nb_sequence == 4
    assert sall.nb_variable == 3

    sall = s1.merge([s3])
    assert sall.nb_sequence == 3
    assert sall.nb_variable == 3


def test_merge_and_Merge(build_seqn):
    s1 = build_seqn
    s2 = build_seqn

    a = s1.merge([s2])
    b = s2.merge([s1])
    v = Merge(s1, s2)

    assert str(a) == str(b)
    assert str(a) == str(v)


def test_imcompatible_merge(build_seqn, build_seq1):
    s1 = build_seqn
    s3 = build_seq1
    try:
        s1.merge(s3)
        assert False
    except:
        assert True


def test_merge_variable(build_seqn, build_seq1):
    s1 = build_seqn
    s2 = build_seqn
    s3 = build_seq1

    sall = s1.merge_variable([s2], 1)  # why 1 ? same result with 2 !
    assert sall.nb_sequence == 2
    assert sall.nb_variable == 6

    # sall =  s1.merge_variable([s3],1)
    # assert sall.nb_sequence == 2
    # assert sall.nb_variable == 3


def test_merge_variable_and_MergeVariable(build_seqn, build_seq1):
    s1 = build_seqn
    s2 = build_seqn
    s3 = build_seq1

    a = s1.merge_variable([s2], 1)
    b = s2.merge_variable([s1], 1)
    v = MergeVariable(s1, s2)

    assert str(a) == str(b)
    assert str(a) == str(v)


def test_cluster_step():
    seq1 = Sequences([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
    assert str(Cluster(seq1, "Step", 1, 2)) == str(seq1.cluster_step(1, 2, True))
    seqn = Sequences([[[1, 2, 3], [1, 3, 1]], [[4, 5, 6], [7, 8, 9]]])
    assert str(Cluster(seqn, "Step", 1, 2)) == str(seqn.cluster_step(1, 2, True))


def test_cluster_limit():
    seq1 = Sequences([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
    assert str(Cluster(seq1, "Limit", 1, [2])) == str(seq1.cluster_limit(1, [2], True))
    seqn = Sequences([[[1, 2, 3], [1, 3, 1]], [[4, 5, 6], [7, 8, 9]]])
    assert str(Cluster(seqn, "Limit", 1, [2, 4, 6])) == str(
        seqn.cluster_limit(1, [2, 4, 6], True)
    )


def test_transcode(build_seqn):
    """This functionality need to be checked.

    See also the vector case!"""
    seq = build_seqn
    assert str(
        seq.transcode(
            1,
            [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0o1, 1, 1, 1, 1, 0, 0],
            False,
        )
    ) == str(
        Transcode(
            seq,
            1,
            [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0o1, 1, 1, 1, 1, 0, 0],
        )
    )


def test_reverse(build_seqn):
    """reverse to be checked. seems to give same output as input"""
    s = build_seqn
    s.reverse()


def test_difference(build_data):
    data = build_data
    assert str(Difference(data, 1)) == str(data.difference(1, False))
    res = Difference(data, 1)
    assert res.cumul_length == 50


def test_cumulate(build_data):
    # see also test_cumulate for more tests
    s = build_data
    res = Cumulate(s)
    assert res.cumul_length == 52

# def test_select_step():
#     """
#     #########################################################################
#     #
#     #  Well-log data; used in Fearnhead and Clifford "On-line Inference for
#     #  Hidden Markov Models via Particle Filters". Measurements of Nuclear-response
#     #  of a well-bore over time. Data from O Ruanaidh, J. J. K. and
#     #  Fitzgerald, W. J. (1996). "Numerical Bayesion Methods Applied to Signal
#     #  Processing". New York: Springer.
#     #
#     #########################################################################
#     """
#     seq1 = Sequences(str(get_shared_data("well_log_filtered.seq")))
#     Plot(seq1, ViewPoint="Data")
#     Plot(seq1)

#     SelectStep(seq1, 1000)
#     Plot(seq1)

#     # Display(seq1, 1, 17, "Gaussian", ViewPoint="SegmentProfile", NbSegmentation=5)
#     Plot(seq1, 1, 17, "Gaussian", ViewPoint="SegmentProfile")

#     # seq20 = Segmentation(seq1, 1, 20, "Gaussian")
#     # seq40 = Segmentation(seq1, 1, 40, "Gaussian")

#     # seq20 = Segmentation(seq1, 1, 20, "Mean")
#     # seq40 = Segmentation(seq1, 1, 40, "Mean")

#     # seq16 = Segmentation(seq1, 1, 16, "Gaussian", NbSegment->"Fixed")

#     vec1 = Vectors(seq1)
#     Plot(vec1)

#     SelectStep(vec1, 1000)
#     Plot(vec1)



if __name__ == "__main__":    

    def build_seqn():
        return Sequences([[[1, 1, 1], [12, 12, 12]], [[2, 2, 2], [22, 23, 24]]])


    def test_shift_seqn(build_seqn):
        s = build_seqn
        shifted = s.shift(1, 2)
        assert shifted[0, 0] == [3, 1, 1]
        assert shifted[0, 1] == [14, 12, 12]
        assert shifted[1, 0] == [4, 2, 2]
        assert shifted[1, 1] == [24, 23, 24]
        shifted1 = s.shift(1, 1)
        assert shifted1[0, 0] == [2, 1, 1]
        assert shifted1[0, 1] == [13, 12, 12]
        assert shifted1[1, 0] == [3, 2, 2]
        assert shifted1[1, 1] == [23, 23, 24]

    test_shift_seqn(build_seqn())