
try:
    from .tools import DISABLE_PLOT, interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import DISABLE_PLOT, interface
    from tools import robust_path as get_shared_data
   
import os

from openalea.sequence_analysis import (
	Sequences,
	SelectVariable,
	Merge,
	Cluster,
	HiddenSemiMarkov,
	Estimate,
	BuildAuxiliaryVariable,
	ExtractData,
	SelectIndividual,
	Plot
)

def test1():
    seq6 = Sequences(str(get_shared_data("pin_laricio_6.seq")))
    seq12 = Sequences(str(get_shared_data("pin_laricio_12.seq")))
    seq18 = Sequences(str(get_shared_data("pin_laricio_18.seq")))
    seq23 = Sequences(str(get_shared_data("pin_laricio_23.seq")))

    seq29 = SelectVariable(Merge(seq6, seq12, seq18, seq23), [1, 2])
    seq30 = Cluster(seq29, "Step", 1, 10)
    # Plot(seq30, 1)

    hsmc29 = HiddenSemiMarkov(str(get_shared_data("pin_laricio_3_gaussian_multivariate.hsc")))
    hsmc30 = Estimate(seq30, "HIDDEN_SEMI-MARKOV", hsmc29)
    # Plot(hsmc30, 1)
    seq31 = BuildAuxiliaryVariable(ExtractData(hsmc30))
    Plot(SelectIndividual(seq31, [95]), ViewPoint="Data")
    
if __name__ == "__main__":
    test1()
