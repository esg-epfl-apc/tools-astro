import subprocess
import sys
import tempfile
from os import makedirs, path
from os.path import join

import numpy as np


class objdict(dict):
    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)


class CRbeam(object):
    default_parameters = dict(
        nparticles=10000,
        z=0.1,
        emax=1e13,
        emin=1e7,
        emin_source=None,  # if None, emin is used
        primary="photon",
        EGMF=0.0,
        lmaxEGMF=5,
        lminEGMF=0.05,
        background=12,
    )

    particles = dict(
        electron=0, positron=1, photon=2, gamma=2, neutron=9, proton=10
    )

    prog_name = "crbeam"

    # Here we always use fixed power law E^{-1} for MC and then adjust weights for the events if needed
    power_law = 1

    def __init__(self, **kwargs):
        self._step = 0
        self._cache = None
        self._size_per_step = 100
        p = objdict(self.default_parameters)
        p.update(kwargs)
        if p.emin_source is None:
            p.emin_source = p.emin

        self._cache_name_prefix = (
            f"{p.primary}_z{p.z}_E_{p.emin:g}_{p.emax:g}_{p.background}"
        )
        if p.emin_source < p.emax:
            self._cache_name_prefix += (
                f"_pow{self.power_law}_Emin{p.emin_source:g}"
            )
        if p.EGMF > 0:
            self._cache_name_prefix += (
                f"_B{p.EGMF:g}_turbL{p.lminEGMF}-{p.lmaxEGMF}"
            )

        self._cache_name_prefix += "_N"

        self._output_dir = self._cache_name_prefix + f"{p.nparticles}"

        p.primary = self.particles[p.primary]

        self.p = p

    @property
    def output_dir(self):
        return self._output_dir

    @property
    def cache(self):
        from cache import CRbeamCache

        if self._cache is None:
            self._cache = CRbeamCache()
        return self._cache

    @property
    def command(self):
        result = f"{self.prog_name} "
        result += " ".join(
            [f"--{key} {value}" for key, value in self.p.items()]
        )
        result += " --output " + self.output_dir
        power = 1 if self.p.emin_source < self.p.emax else -1000
        result += f" --power {power}"
        if self.p.EGMF > 0:
            # use turbulent magnetic field with unique random orientation for all particles
            result += " -mft -mfr"
        return result

    @property
    def output_path(self):
        return join(self.output_dir, "z0")

    def run(self, force_overwrite):
        command_suffix = ""
        if force_overwrite:
            command_suffix = " --overwrite"
        if path.isdir(self.output_path) and not force_overwrite:
            return
        logs_dir = "logs"
        makedirs(logs_dir, exist_ok=True)
        log_file = f"{logs_dir}/{self.output_dir}.log"
        with tempfile.NamedTemporaryFile() as script_file:
            with open(script_file.name, "wt") as out_script:
                print("#!/bin/bash", file=out_script)
                print(
                    self.command + command_suffix,
                    "2>&1 1>" + log_file,
                    file=out_script,
                )
            subprocess.run(["bash", script_file.name])
        return log_file

    def run_cached(self, overwrite_local_cache=False):
        import shutil

        output_root = self.output_dir
        if path.isdir(output_root):
            if not overwrite_local_cache:
                return self.output_path
            shutil.rmtree(output_root)

        size = self.cache.get_cache_size(self._cache_name_prefix)
        print("s3 data size: ", size)
        skip_paths = []
        makedirs(self.output_path, exist_ok=False)
        if size < self.p.nparticles:
            try:
                self.p.nparticles = self.p.nparticles - size
                self.run(force_overwrite=True)
                skip_paths.append(self.upload_cache())
            finally:
                self.p.nparticles = self.p.nparticles + size
        else:
            self.p.nparticles = size

        if size > 0:  # append results from s3 cache to the local cach file
            self.cache.load_results(
                self.output_path + "/photon",
                self._cache_name_prefix,
                skip_paths=skip_paths,
            )

        return self.output_path

    def start_multistage_run(self, overwrite_local_cache=False, n_steps=10):
        """
        Initialize environment for running different parts of simulation in subsiquent calls of simulation_step
        This function is useful if we show progress by running code in different jupyter cells
        :param overwrite_local_cache: remove local cache
        :param n_steps: number intermediate steps
        :return: True if further simulation is needed, False if data is already available
        """
        import shutil

        output_root = self.output_dir
        if path.isdir(output_root):
            if not overwrite_local_cache:
                return False
            shutil.rmtree(output_root)

        size = self.cache.get_cache_size(self._cache_name_prefix)
        print("s3 data size: ", size)
        if size >= self.p.nparticles:
            makedirs(self.output_path, exist_ok=False)
            self.cache.load_results(
                self.output_path + "/photon", self._cache_name_prefix
            )
            self.p.nparticles = size
            return False
        self.p.nparticles = self.p.nparticles - size
        self._size_per_step = max(
            100, int(np.ceil((self.p.nparticles) / n_steps))
        )
        self._step = 0
        return True

    def simulation_step(self):
        size = self._size_per_step * self._step
        target_nparticles = self.p.nparticles
        if size >= target_nparticles:
            return False
        save_output_dir = self._output_dir
        try:
            adjusted_size_per_step = min(
                self._size_per_step, target_nparticles - size
            )
            self.p.nparticles = adjusted_size_per_step
            self._step += 1
            self._output_dir = f"{self._output_dir}_step{self._step}"
            self.run(force_overwrite=True)
        finally:
            self.p.nparticles = target_nparticles
            self._output_dir = save_output_dir

        return size + adjusted_size_per_step < target_nparticles

    def end_multistep_run(self):
        makedirs(self.output_path, exist_ok=False)
        with tempfile.NamedTemporaryFile() as script_file:
            with open(script_file.name, "wt") as task:
                print(
                    f"cat {self.output_dir}_step*/z0/photon >"
                    f" {self.output_path}/photon",
                    file=task,
                )
            subprocess.run(["bash", script_file.name])
        try:
            skip_paths = [self.upload_cache()]
            self.cache.load_results(
                self.output_path + "/photon",
                self._cache_name_prefix,
                skip_paths=skip_paths,
            )
            self.p.nparticles = self.cache.get_cache_size(
                self._cache_name_prefix
            )
        except Exception as ex:
            print(ex, file=sys.stderr)
        return self.output_path

    def upload_cache(self):
        source_file = self.output_path + "/photon"
        if not path.isfile(source_file):
            if path.isdir(self.output_path):
                print("empty result")
                return None
            else:
                assert False, "result dir not found"
        obj_name = self._cache_name_prefix + f"{self.p.nparticles}"
        print("saving", source_file, "to s3 as", obj_name)
        return self.cache.save(obj_name, source_file)

    def remove_cache(self):
        from cache import CRbeamCache

        c = CRbeamCache()
        prefix = (
            self._cache_name_prefix
        )  # todo: use metadata when metadata reading is implemented
        c.detete_results(prefix)
