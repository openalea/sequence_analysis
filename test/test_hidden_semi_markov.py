# -*- coding: utf-8 -*-
"""hidden semi markov data structure  tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id$"


#import openalea.stat
from openalea.stat_tool import _stat_tool, set_seed
from openalea.sequence_analysis import (_sequence_analysis, 
                                        Estimate,
                                        seq_map)

from openalea.sequence_analysis.hidden_semi_markov import HiddenSemiMarkov
from openalea.sequence_analysis.simulate import Simulate
from openalea.sequence_analysis.data_transform import Thresholding

from openalea.stat_tool.data_transform import *
from openalea.stat_tool.cluster import Transcode, Cluster

from openalea.stat_tool import Display, Save
try:
    from .tools import DISABLE_PLOT, interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import DISABLE_PLOT, interface
    from tools import robust_path as get_shared_data

import os

import pytest 
from dataclasses import dataclass
from typing import Any


@dataclass
class HSM:
    data: Any
    hsm: Any

@pytest.fixture
def HSMData():
     data = HiddenSemiMarkov(get_shared_data( "wij1.hsc"))
     hsm = HiddenSemiMarkov(get_shared_data('test_hidden_semi_markov.dat'))
     return HSM(data, hsm)

def _remove_file(filename):
    """alias to remove a file"""
    try:
        os.remove(filename)
    except:
        pass

class TestHiddenSemiMarkov():
    """
    Test HiddenSemiMarkov class
    """
    def init(self, HSMData):        
        data = HSMData.data
        assert data
        self.data = data
        hsm = HSMData.hsm
        assert hsm
        self.hsm = hsm
        self.structure = HiddenSemiMarkov



    def constructor_from_file(self):
        """Test constructor from file"""
        if self.filename == None:
            return None
        else:
            c = self.structure(self.filename)
            assert c
            return c

    def constructor_from_file_failure(self):
        """run constructor with filename argument"""
        try:
            _h = self.structure("whatever_wrong_filename.txt")
            assert False
        except Exception:
            assert True

    def print_data(self):
        """test that print command exists"""
        print(self.data)

    def display(self):
        """check that .display/Display is callable"""
        data = self.data
        data.display()
        Display(data)
        assert data.display() == Display(data)

    def display_versus_ascii_write(self):
        """check that display is equivalent to ascii_write"""
        assert Display(self.data) == self.data.ascii_write(False)

    def display_versus_str(self):
        """check that display and str are equivalent"""
        data = self.data
        s = str(data)
        assert Display(data) == s

    def plot(self):
        """run plotting routines"""
        # if DISABLE_PLOT == False:
        self.data.plot()

    def save(self, Format=None, skip_reading=False):
        """In the Vector case, Format should be set to Data.
        :param skip_reading: some class do not have Filename Constructor;
            skip_reading can be set to False to prevent code to be run.

        .. todo:: This is surely a bug. to be checked"""

        c1 = self.data
        # _remove_file('test1.dat')
        # _remove_file('test2.dat')

        if Format is None:
            c1.save("test1.dat")
            Save(c1, "test2.dat")
        else:
            c1.save("test1.dat", Format="Data")
            Save(c1, "test2.dat", Format="Data")

        if skip_reading:
            pass
        else:
            c1_read = self.structure("test1.dat")
            c2_read = self.structure("test2.dat")

            assert c1 and c1_read and c2_read
            assert str(c1_read) == str(c2_read)

        _remove_file("test1.dat")
        _remove_file("test2.dat")

    def plot_write(self):
        h = self.data
        h.plot_write("test", "title")
        _remove_file("test.print")
        _remove_file("test.plot")
        _remove_file("test0.dat")

    def file_ascii_write(self):
        h = self.data
        h.file_ascii_write("test.dat", True)
        _remove_file("test.dat")

    def file_ascii_data_write(self):
        h = self.data
        h.file_ascii_data_write("test.dat", True)
        _remove_file("test.dat")

    def spreadsheet_write(self):
        h = self.data
        h.spreadsheet_write("test.dat")
        _remove_file("test.dat")

    def survival_ascii_write(self):
        d = self.data
        d.survival_ascii_write()

    def survival_plot_write(self):
        d = self.data
        d.survival_plot_write("test", "test")

    def survival_file_ascii_write(self):
        d = self.data
        d.survival_file_ascii_write()

    def survival_spreadsheet_write(self):
        d = self.data
        d.survival_spreadsheet_write("test.xsl")
        _remove_file("test.xsl")

    def simulate(self, N=-1):
        """Test the simulate method"""
        set_seed(0)
        if N == -1:
            N = self.N
        m = self.data
        s = m.simulate(N)
        s2 = Simulate(m, N)
        assert len(s) == N
        assert len(s2) == N
        assert str(s2)
        assert str(s)
        return s

    def test_len(self):
        raise NotImplementedError()

    def test_extract(self):
        raise NotImplementedError()

    def test_extract_data(self):
        raise NotImplementedError()

    def test_constructor_from_file2(self):
        # self.init()
        hmc = HiddenSemiMarkov(get_shared_data("test_hidden_markov.hmc"))
        assert hmc

    def test_constructor_from_file_nonparametric1(self):
        """Read HSM model from a file with 1st nonparametric variable
        with observation distribution ending by 1e-05 """
        hmc = HiddenSemiMarkov(get_shared_data("test_hidden_markov_non-parametric1.hmc"))
        assert hmc

    def test_constructor_from_file_failure(self, HSMData):
        self.init(HSMData)
        self.constructor_from_file_failure()

    def test_print(self, HSMData):
        self.init(HSMData)
        self.print_data()

    def test_len(self, HSMData):
        """
        test_len Not applicable
        Implemented since deriving from interface
        """
        self.init(HSMData)
        pass

    def test_display(self, HSMData):
        self.init(HSMData)
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_plot(self, HSMData):
        self.init(HSMData)
        self.plot()

    def test_save(self, HSMData):
        self.init(HSMData)
        self.save(skip_reading=True)

    def test_plot_write(self, HSMData):
        self.init(HSMData)
        self.plot_write()

    def test_file_ascii_write(self, HSMData):
        self.init(HSMData)
        self.file_ascii_write()

    def test_spreadsheet_write(self, HSMData):
        self.init(HSMData)
        self.spreadsheet_write()

    def simulate(self, HSMData):
        self.init(HSMData)
        sm = self.hsm
        assert sm.simulation_nb_elements(1, 10000, True)
        s = Simulate(sm,1, 10000, True)
        return s
    
    def test_simulate(self, HSMData):
        s = self.simulate(HSMData)
        assert s
        
    def test_thresholding(self, HSMData):
        self.init(HSMData)
        a = self.data.thresholding(0.01)
        b = Thresholding(self.data, MinProbability=0.01)
        assert str(a)==str(b)

    def extract_data(self, HSMData):
        self.init(HSMData)
        assert self.hsm.extract_data() is None
        nb_seq = 10
        seq_length = 50
        set_seed(0)
        seq = self.hsm.simulation_nb_sequences(nb_seq, seq_length, True)
        # Discard state
        obs = seq.select_variable([1], False)
        hsm_estim = Estimate(obs, "HIDDEN_SEMI-MARKOV", self.hsm, NbIteration=3)
        return hsm_estim

    def test_extract_data(self, HSMData):
        hsm_estim = self.extract_data(HSMData)
        assert hsm_estim.extract_data()

    def test_nb_output_process(self, HSMData):
        """Test consistency of the number of output processes"""        
        self.init(HSMData)
        s = self.simulate(HSMData)
        assert self.hsm.nb_output_process == (s.nb_variable - 1)

    def test_ascii_write(self, HSMData):
        """Test consistency of text representation"""
        self.init(HSMData)
        assert str(self.data) == self.data.ascii_write(False) 

    def test_simulation_histogram(self, HSMData):
        """Test simulation from a distribution of sequence lengths"""
        self.init(HSMData)
        s1 = self.simulate(HSMData)
        s2 = self.hsm.simulation_histogram(s1.extract_length(), False, False)
        assert str(s1.extract_length()) == str(s2.extract_length())

    def test_state_sequence_computation(self, HSMData):
        """Test state sequence restoration"""  
        h = self.extract_data(HSMData)
        s = h.extract_data()
        r = r = h.state_sequence_computation(s.select_variable([1], False), True)
        assert self.hsm.nb_output_process == r.nb_variable-1

    def test_extract(self, HSMData):
        """Test state sequence restoration"""  
        h = self.extract_data(HSMData)
        for k in range(h.get_nb_state()):
            if k < h.get_nb_state()-1:
                assert h.extract(seq_map['Sojourn'],k,0)
            v = 1 # Output process 1
            assert v <= h.nb_output_process 
            dist = h.extract(seq_map['Observation'],v,k)
            assert dist 
            m = dist.get_alloc_nb_value-1
            for o in range(m):
                assert h.extract(seq_map['Sojourn'],v,o)
                assert h.extract(seq_map['FirstOccurrence'],v,o)
                assert h.extract(seq_map['Recurrence'],v,o)
                assert h.extract(seq_map['NbRun'],v,o)
                assert h.extract(seq_map['NbOccurrence'],v,o)
"""
TODO: 
hsm.divergence_computation
hsm.nb_iterator
"""


if __name__ == "__main__":
    
    class HSM:
        def __init__(self, data, hsm):
            self.data = data
            self.hsm = hsm
    
    def HSMData():
        data = HiddenSemiMarkov(get_shared_data( "wij1.hsc"))
        hsm = HiddenSemiMarkov(get_shared_data('test_hidden_semi_markov.dat'))
        return HSM(data, hsm)

    T = TestHiddenSemiMarkov()
    T.init(HSMData())
    T.test_constructor_from_file_failure(HSMData())
    T.test_constructor_from_file2()
    T.test_constructor_from_file_nonparametric1()
    T.test_print(HSMData())
    T.test_display(HSMData())
    T.test_len(HSMData())
    T.test_plot(HSMData())
    T.test_save(HSMData())
    T.test_plot_write(HSMData())
    T.test_file_ascii_write(HSMData())
    T.test_spreadsheet_write(HSMData())
    T.test_simulate(HSMData())
    T.test_thresholding(HSMData())
    T.test_extract(HSMData())
    T.test_extract_data(HSMData())
    # T.test_extract_histogram(HSMData())
    # T.test_state_sequence_computation(HSMData())