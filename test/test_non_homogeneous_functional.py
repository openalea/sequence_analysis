#########################################################################
#
#  Floraison des tiges de vanille decrites du sommet vers la base.
#  deux echantillons : apex mort (m) ou decapite (d).
#
#  VARIABLE 1 : non-fleuri (0) / fleuri (1).
#
#########################################################################
from openalea.sequence_analysis import Sequences, Estimate, NonhomogeneousMarkov
from openalea.sequence_analysis import ComputeSelfTransition, Plot

try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data
    

# mc_m = Estimate(seq_m, "NONHOMOGENEOUS_MARKOV", "MONOMOLECULAR", "VOID")
# Plot(mc_m, "SelfTransition")
# Plot(mc_m, "Intensity")

