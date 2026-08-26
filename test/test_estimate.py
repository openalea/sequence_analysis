"""tests on the method estimate

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: finalise
"""

__revision__ = "$Id$"

import pytest

from openalea.stat_tool.vectors import Vectors
from openalea.stat_tool.data_transform import ExtractHistogram
from openalea.stat_tool import set_seed

try:
    from .tools import DISABLE_PLOT, interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import DISABLE_PLOT, interface
    from tools import robust_path as get_shared_data

from openalea.sequence_analysis import (
    Sequences,
    Merge,
    RemoveRun,
    SegmentationExtract,
    LengthSelect,
    Estimate,
    _SemiMarkovData,
    HiddenVariableOrderMarkov,
    HiddenSemiMarkov
)

_seq1 = Sequences(get_shared_data("dupreziana_20a2.seq"))

seq2 = RemoveRun(_seq1, 1, 0, "End")
seq3 = Sequences(get_shared_data("dupreziana_40a2.seq"))
seq4_0 = RemoveRun(seq3, 2, 0, "End")
seq4 = SegmentationExtract(seq4_0, 1, 2)
seq5 = Sequences(get_shared_data("dupreziana_60a2.seq"))
seq6_0 = RemoveRun(seq5, 2, 0, "End")
seq6 = LengthSelect(SegmentationExtract(seq6_0, 1, 2), 1, Mode="Reject")
seq7 = Sequences(get_shared_data("dupreziana_80a2.seq"))
seq8_0 = RemoveRun(seq7, 2, 0, "End")
seq8 = SegmentationExtract(seq8_0, 1, 2)
seq10 = Merge(seq2, seq4, seq6, seq8)


@pytest.fixture
def create_data_estimate_histogram():
    seq0 = Sequences(get_shared_data("chene_sessile_15pa.seq"))
    return Vectors(seq0)


def test_estimate_mixture(create_data_estimate_histogram):
    set_seed(0)
    mixt20 = Estimate(
        ExtractHistogram(create_data_estimate_histogram, 2),
        "MIXTURE",
        "NB",
        "NB",
        "NB",
        "NB",
        NbComponent="Estimated",
    )
    assert mixt20.nb_component == 2


def test_estimate_mixture2(create_data_estimate_histogram):
    set_seed(0)
    mixt20 = Estimate(
        ExtractHistogram(create_data_estimate_histogram, 5),
        "MIXTURE",
        "NB",
        "NB",
        "NB",
        "NB",
        NbComponent="Estimated",
    )
    assert mixt20.nb_component == 4


sequence = seq10
estimate_type = "VARIABLE_ORDER_MARKOV"


def test_estimate():
    set_seed(0)
    mc10 = Estimate(
        sequence,
        estimate_type,
        "Ordinary",
        MaxOrder=5,
        GlobalInitialTransition=True,
    )


def test_estimate1():
    set_seed(0)
    mc11 = Estimate(
        sequence,
        estimate_type,
        "Ordinary",
        MaxOrder=5,
        GlobalInitialTransition=False,
    )


def test_estimate2():
    set_seed(0)
    mc12 = Estimate(
        sequence,
        estimate_type,
        "Ordinary",
        Algorithm="LocalBIC",
        Threshold=10.0,
        MaxOrder=5,
        GlobalInitialTransition=False,
        GlobalSample=False,
    )


def test_estimate3():
    set_seed(0)
    mc13 = Estimate(
        sequence,
        estimate_type,
        "Ordinary",
        Algorithm="Context",
        Threshold=1.0,
        MaxOrder=5,
        GlobalInitialTransition=False,
        GlobalSample=False,
    )


def test_estimate4():
    set_seed(0)
    for Algorithm in ["CTM_BIC", "CTM_KT", "Context"]:
        mc13 = Estimate(
            sequence,
            estimate_type,
            "Ordinary",
            Algorithm=Algorithm,
            MaxOrder=5,
            GlobalInitialTransition=False,
            GlobalSample=False,
        )


def test_estimate_error1():
    """test that Estimator and Algorith=CTM_KT are incompatible"""
    set_seed(0)
    try:
        mc13 = Estimate(
            sequence,
            estimate_type,
            "Ordinary",
            Algorithm="CTM_KT",
            Estimator="Laplace",
            MaxOrder=5,
            GlobalInitialTransition=False,
            GlobalSample=False,
        )
        assert False
    except:
        assert True


def Test_Estimate_VARIABLE_ORDER_MARKOV_from_markovian():
    set_seed(0)
    mc11 = Estimate(
        seq10,
        "VARIABLE_ORDER_MARKOV",
        "Ordinary",
        MaxOrder=5,
        GlobalInitialTransition=False,
    )
    mc2 = Estimate(seq2, "VARIABLE_ORDER_MARKOV", mc11, GlobalInitialTransition=False)


@pytest.fixture
def create_sequence_estimate_variable_order_markov():
    return Sequences(get_shared_data("sequences1.seq"))

@pytest.fixture
def create_hvom_estimate_hidden_variable_order_markov():
    return HiddenVariableOrderMarkov(get_shared_data("dupreziana21.hc"))

"""
def test_estimate_hidden_variable_order_markov(
    create_sequence_estimate_variable_order_markov,
    create_hvom_estimate_hidden_variable_order_markov,
):
    set_seed(0)
    hmc_estimated = Estimate(
        create_sequence_estimate_variable_order_markov,
        "HIDDEN_VARIABLE_ORDER_MARKOV",
        create_hvom_estimate_hidden_variable_order_markov,
        GlobalInitialTransition=True,
        NbIteration=80,
    )
    assert hmc_estimated
"""

@pytest.fixture
def create_data_estimate_hidden_semi_markov():
    return HiddenSemiMarkov(get_shared_data("wij1.hsc"))


@pytest.fixture
def create_sequence_estimate_hidden_semi_markov():
    return Sequences(get_shared_data("wij1.seq"))

@pytest.fixture
def create_sequence_simulate_hidden_semi_markov():
    hsm = HiddenSemiMarkov(get_shared_data('test_hidden_semi_markov.dat'))
    set_seed(0)
    nb_seq = 30
    seq_length = 100
    seq = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
    return seq


def test_estimate_hidden_semi_markov(
    create_data_estimate_hidden_semi_markov, create_sequence_estimate_hidden_semi_markov
):
    """Estimate hidden semi-Markov model from initial number of states"""
    set_seed(0)
    # data is a hsm class
    Estimate(
        create_sequence_estimate_hidden_semi_markov,
        "HIDDEN_SEMI-MARKOV",
        # create_data_estimate_hidden_semi_markov,
        "Ordinary", 3, "LeftRight", NbIteration=300  
    )

def test_estimate_init_model_hidden_semi_markov(
    create_data_estimate_hidden_semi_markov, create_sequence_estimate_hidden_semi_markov
):
    """Estimate hidden semi-Markov model from initial model read from file"""
    set_seed(0)
    # data is a hsm class
    Estimate(
        create_sequence_estimate_hidden_semi_markov,
        "HIDDEN_SEMI-MARKOV",
        create_data_estimate_hidden_semi_markov, 
        NbIteration=300  
    )

def test_estimate_hidden_semi_markov_wrong_nb_variable(
    create_data_estimate_hidden_semi_markov, create_sequence_simulate_hidden_semi_markov
):
    """Estimate HiddenSemiMarkov with initial model and wrong number of variables"""
    set_seed(0)
    # data is a hsm class
    create_sequence_simulate_hidden_semi_markov.set_type_to_int(1);
    assert(create_sequence_simulate_hidden_semi_markov.get_type(1) == 0)
    try:
        Estimate(
            create_sequence_simulate_hidden_semi_markov,
            "HIDDEN_SEMI-MARKOV",
            create_data_estimate_hidden_semi_markov,
            NbIteration=300  
        )
    except:
        assert True
    else:
        assert False

def test_estimate_hidden_semi_markov_wrong_type(
    create_sequence_simulate_hidden_semi_markov
):
    """Estimate HiddenSemiMarkov with initial model and wrong variable type"""
    set_seed(0)
    # data is a hsm class
    try:
        Estimate(
        create_sequence_simulate_hidden_semi_markov,
        "HIDDEN_SEMI-MARKOV",
        "Ordinary", 3, "LeftRight", NbIteration=300  
        )
    except:
        assert True
    else:
        assert False
"""
def test_estimate_semi_markov():
    set_seed(0)
    sequence = Sequences(get_shared_data("wij1.seq"))
    Estimate(sequence, "SEMI-MARKOV", "Ordinary")
"""

def test_estimate_time_events():
    """test not yet implemented"""
    pass




if __name__ == "__main__":
    def create_data_estimate_histogram():
        seq0 = Sequences(get_shared_data("chene_sessile_15pa.seq"))
        return Vectors(seq0)

    def create_data_estimate_hidden_semi_markov():
        return HiddenSemiMarkov(get_shared_data("wij1.hsc"))

    def create_sequence_estimate_variable_order_markov():
        return Sequences(get_shared_data("sequences1.seq"))

    def create_hvom_estimate_hidden_variable_order_markov():
        return HiddenVariableOrderMarkov(get_shared_data("dupreziana21.hc"))

    def create_sequence_estimate_hidden_semi_markov():
        return Sequences(get_shared_data("wij1.seq"))

    def create_sequence_simulate_hidden_semi_markov():
        hsm = HiddenSemiMarkov(get_shared_data('test_hidden_semi_markov.dat'))
        set_seed(0)
        nb_seq = 30
        seq_length = 100
        seq = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
        return seq

    def test_estimate_mixture():
        set_seed(0)
        mixt20 = Estimate(
            ExtractHistogram(create_data_estimate_histogram(), 2),
            "MIXTURE",
            "NB",
            "NB",
            "NB",
            "NB",
            NbComponent="Estimated",
        )
        assert mixt20.nb_component == 2
    
    def test_estimate_hidden_semi_markov():
        set_seed(0)
        # data is a hsm class
        Estimate(
            create_sequence_estimate_hidden_semi_markov(),
            "HIDDEN_SEMI-MARKOV",
            create_data_estimate_hidden_semi_markov(),
        )
        # hsmc_est = Estimate(obs, "HIDDEN_SEMI-MARKOV", hsm, Nbiteration=300) 
    
    def test_estimate_hidden_semi_markov_wrong_type():
        """Estimate HiddenSemiMarkov with initial model and wrong variable type"""
        set_seed(0)
        # data is a hsm class
        try:
            Estimate(
            create_sequence_simulate_hidden_semi_markov(),
            "HIDDEN_SEMI-MARKOV",
            "Ordinary", 3, "LeftRight", NbIteration=300  
        )
        except:
            assert True
        else:
            assert False

    def test_estimate_hidden_semi_markov_wrong_nb_variable():
        """Estimate HiddenSemiMarkov with initial model and wrong number of variables"""
        set_seed(0)
        # data is a hsm class
        seq = create_sequence_simulate_hidden_semi_markov()
        seq.set_type_to_int(1);
        assert(seq.get_type(1) == 0)
        try:
            Estimate(
                seq,
                "HIDDEN_SEMI-MARKOV",
                create_data_estimate_hidden_semi_markov(),
                NbIteration=300  
            )
        except:
            assert True
        else:
            assert False
    def test_estimate_mixture2():
        set_seed(0)
        mixt20 = Estimate(
            ExtractHistogram(create_data_estimate_histogram(), 5),
            "MIXTURE",
            "NB",
            "NB",
            "NB",
            "NB",
            NbComponent="Estimated",
        )
        assert mixt20.nb_component == 4


    
    """
    test_estimate_mixture()
    test_estimate_mixture2()
    test_estimate()
    test_estimate1()
    test_estimate2()
    test_estimate3()
    test_estimate4()
    test_estimate_error1()"""
    # Test_Estimate_VARIABLE_ORDER_MARKOV_from_markovian()
    test_estimate_hidden_semi_markov()
    test_estimate_hidden_semi_markov_wrong_nb_variable()
    """
    test_estimate_hidden_semi_markov()
    test_estimate_semi_markov()
    """
