import sscha.Cluster as Cluster
import sys, os

"""
Define a local cluster class.
This allows to mock the cluster class and run the code locally, but by
using the same interface as the cluster class and a job scheduler like SLURM.
"""


class LocalCluster(Cluster.Cluster):
    def ExecuteCMD(self, cmd, *args, on_cluster = False, **kwargs):
        """
        Execute a command in the local machine.
        """

        # Override the value of on_cluster
        return super().ExecuteCMD(cmd, *args, on_cluster = False, **kwargs)

    def copy_file(self, source, destination, server_source = False, server_dest = False, **kwargs):
        """
        Copy the files ignoring if the cluster is used.
        """

        return super().copy_file(source, destination, server_source = False, server_dest = False, **kwargs)
