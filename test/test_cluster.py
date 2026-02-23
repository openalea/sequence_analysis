""" Cluster tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: check the AddVariable option (sequences) and sequences cases
"""
__revision__ = "$Id$"

import os
from openalea.sequence_analysis.sequences import Sequences
#from openalea.sequence_analysis.semi_markov import SemiMarkov
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.vectors import Vectors

try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data


class _Cluster():
    """Test class to test cluster function and classes

    create_data, cluster_step and cluster_limit funciton will be required
    """

    def init(self):
        self.data = None

    def create_data(self):
        self.init()
        raise NotImplemented

    def test_cluster_step(self):
        raise NotImplemented

    def test_cluster_limit(self):
        raise NotImplemented

class _HistoCase(_Cluster):
    """
    inherits from _cluster and implements the cluster_limit and cluster_step
    functions.

    In addition, classes that inherits from _HistoCase must implement
    cluster_information
    """
    def init(self):
        super().init()
        self.data = None

    def test_cluster_step(self):
        self.create_data()
        if not(self.data is None):
            cluster1 = Cluster(self.data, "Step", 2)
            cluster2 = self.data.cluster_step(2)
            assert str(cluster1) == str(cluster2)

    def test_cluster_limit(self):
        self.create_data()
        if not(self.data is None):
            cluster1 = Cluster(self.data, "Limit", [2, 4, 6, 8, 10])
            cluster2 = self.data.cluster_limit([2, 4, 6, 8, 10])
            assert str(cluster1) == str(cluster2)

    def test_cluster_information(self):
        self.create_data()
        if not(self.data is None):
            cluster1 = Cluster(self.data, "Information", 0.8)
            cluster2 = self.data.cluster_information(0.8)
            assert str(cluster1) == str(cluster2)


class TestHistogram(_HistoCase):

    def init(self):
        super().init()

    def create_data(self):
        self.init()
        return Histogram(get_shared_data("data/fagus1.his"))


class TestConvolution( _HistoCase):

    def init(self):
        super().init()

    def create_data(self):
        self.init()
        conv = Convolution(get_shared_data("data/test_convolution1.conv"))
        return conv.simulate(1000)


class TestCompound(_HistoCase):

    def init(self):
        super().init()

    def create_data(self):
        self.init()
        comp = Compound(get_shared_data("data/test_compound1.cd"))
        return comp.simulate(1000)


class TestVectorsn(_Cluster):

    def init(self):
        super().init()

    def create_data(self):
        self.init()
        data = Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        self.data = data
        return self.data

    def test_cluster_step(self):
        self.create_data()
        cluster1 = self.data.cluster_step(1, 2)
        cluster2 = Cluster(self.data, "Step", 1, 2)
        assert str(cluster1) == str(cluster2)

    def test_cluster_limit(self):
        self.create_data()
        cluster1 = self.data.cluster_limit(1, [2, 4, 6])
        cluster2 = Cluster(self.data, "Limit", 1, [2, 4, 6])
        assert str(cluster1) == str(cluster2)

class TestVectors1(_Cluster):

    def init(self):
        super().init()

    def create_data(self):
        self.init()
        data = Vectors([[1], [1], [4]])
        self.data = data
        return self.data

    def test_cluster_step(self):
        self.create_data()
        cluster1 = self.data.cluster_step(1, 2)
        cluster2 = Cluster(self.data, "Step", 2)
        assert str(cluster1) == str(cluster2)

    def test_cluster_limit(self):
        self.create_data()
        cluster1 = self.data.cluster_limit(1, [2, 4, 6])
        cluster2 = Cluster(self.data, "Limit",  [2, 4, 6])
        assert str(cluster1) == str(cluster2)


class TestSequences1(_Cluster):

    def init(self):
        super().init()

    def create_data(self):
        self.init()
        data = Sequences(get_shared_data("data/sequences1.seq"))
        self.data = data
        return data

    def test_cluster_step(self):
        self.create_data()
        mode = False
        cluster1 = self.data.cluster_step(1, 2, mode)
        cluster2 = Cluster(self.data, "Step", 2)
        assert str(cluster1) == str(cluster2)

    def test_cluster_limit(self):
        self.create_data()
        print(self.data.nb_variable)
        cluster1 = self.data.cluster_limit(1,[2], False)
        cluster2 = Cluster(self.data,"Limit", [2] , AddVariable=False)
        assert str(cluster1) == str(cluster2)


class TestSequencesn(_Cluster):

    def init(self):
        super().init()

    def create_data(self):
        self.init()
        data = Sequences(get_shared_data("data/sequences2.seq"))
        self.data = data
        return data

    def test_cluster_step(self):
        self.create_data()
        mode = True
        cluster1 = self.data.cluster_step(1, 2, mode)
        cluster2 = Cluster(self.data, "Step", 1, 2)
        assert str(cluster1) == str(cluster2)

    def test_cluster_limit(self):
        self.create_data()
        cluster1 = self.data.cluster_limit(1, [2 ], True)
        cluster2 = Cluster(self.data, "Limit", 1, [2])
        assert str(cluster1) == str(cluster2)

if __name__ ==  "__main__":
    TV1 = TestVectors1()
    TV1.test_cluster_limit()
    TV1.test_cluster_step()
    TVn = TestVectorsn()
    TVn.test_cluster_limit()
    TVn.test_cluster_step()
    TS1 = TestSequences1()
    TS1.test_cluster_limit()
    TS1.test_cluster_step()
    TSn = TestSequencesn()
    TSn.test_cluster_limit()
    TSn.test_cluster_step()
    TC = TestConvolution()
    TC.test_cluster_information()
    TC.test_cluster_limit()
    TC.test_cluster_step()
    TC2 = TestCompound()
    TC2.test_cluster_information()
    TC2.test_cluster_limit()
    TC2.test_cluster_step()
    H = TestHistogram()
    H.test_cluster_information()
    H.test_cluster_limit()
    H.test_cluster_step()
    