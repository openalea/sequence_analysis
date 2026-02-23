"""tests on the method AddAbsorbingRun

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr 

.. todo:: markov case ? 
"""
__revision__ = "$Id$"

from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.semi_markov import SemiMarkov
from openalea.sequence_analysis.data_transform import AddAbsorbingRun
from openalea.stat_tool.distribution import set_seed

try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data


class _AddAbsorbingRun():
    """
    a main class to test the AddAbsorbingrun function on different type of 
    data structure.
    
    """

    def init(self):
        self.data = None 
        self.max_length = -1
        self.MAX_RUN_LENGTH = 20 # hardcoded values in CPP code
        
    def create_data(self):
        self.init()
        self.data = None
        return self.data
        
    def test_max_length(self):
        self.create_data()
        seq = self.data
        assert seq.max_length == self.max_length
        
    def test_boost_versus_module(self):
        self.create_data()
        seq = self.data
        sequence_length = -1
        run_length = 6
        
        boost = seq.add_absorbing_run(sequence_length, seq.max_length+run_length)
        module1 = AddAbsorbingRun(seq, 
                                  SequenceLength=sequence_length,
                                  RunLength=seq.max_length+run_length)
        module2 = AddAbsorbingRun(seq, RunLength=seq.max_length+run_length)

        assert str(module1)==str(boost)
        assert str(module2)==str(boost)
        
    def test_no_arguments(self):
        self.create_data()
        seq = self.data
        assert AddAbsorbingRun(seq)
    
    def test_wrong_run_length(self):
        self.create_data()
        seq = self.data       
        try:
            #second arguments must be less than MAX_RU_LENGTH
            _res = AddAbsorbingRun(seq, -1, self.MAX_RUN_LENGTH + 1)
            assert False
        except Exception:
            assert True

    def test_wrong_sequence_length(self):
        self.create_data()
        seq = self.data       
        try:
            #second arguments must be less than MAX_RU_LENGTH
            _res = AddAbsorbingRun(seq, self.max_length -1, -1)
            assert False
        except Exception:
            assert True

class Test_AddAbsorbingRun_Sequences(_AddAbsorbingRun):
    """sequences case"""
    def init(self):
        _AddAbsorbingRun.init(self)
        self.max_length = 30

    def create_data(self):
        self.init()
        seq = Sequences(str(get_shared_data('data/sequences1.seq')))
        self.data = seq
        return self.data
    
class Test_AddAbsorbingRun_SemiMarkov(_AddAbsorbingRun):
    """semi markov case"""
    def init(self):
        _AddAbsorbingRun.init(self)
        self.max_length = 1000

    def create_data(self):
        self.init()
        markov = SemiMarkov(str(get_shared_data('data/test_semi_markov.dat')))
        semi_markov_data = markov.simulation_nb_elements(1, 1000, True)
        self.data = semi_markov_data
        return self.data
    
if __name__ == "__main__":
    T1 = Test_AddAbsorbingRun_Sequences()
    T1.test_max_length()
    T1.test_boost_versus_module()
    T1.test_no_arguments()
    T1.test_wrong_run_length()
    T1.test_wrong_sequence_length()
    T2 = Test_AddAbsorbingRun_SemiMarkov()
    T2.test_max_length()
    T2.test_boost_versus_module()
    T2.test_no_arguments()
    T2.test_wrong_run_length()
    T2.test_wrong_sequence_length()