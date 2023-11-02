# quick-and-dirty update to DiameterIteration class from opencgs to adjust it for VGF
# this is just a copy of the original code with a few modification
import yaml
from copy import deepcopy
import numpy as np
import shutil
import matplotlib.pyplot as plt
from opencgs.sim import SteadyStateSim, Simulation


class DiameterIterationVGF(Simulation):
    """In this class functionality for simulations with iterative
    crystal diameter computation is collected."""
    def __init__(
        self,
        geo,
        geo_config,
        sim,
        sim_config,
        mat_config,
        config_update={},
        T_tp=505,
        h_min=0.001,
        h_max=0.02,
        max_iterations=10,
        dT_max=0.01,
        sim_name="opencgs-diameter-iteration",
        base_dir="./simdata",
        with_date=True,
        metadata="",
        **_,
    ):
        """Run a simulation with iterative crystal diameter computation,
        consisting of multiple steady-state simulations.

        Args:
            geo (function): function for geometry generation
            geo_config (dict): configuration for geo
            sim (function): function for simulation setup
            sim_config (dict): configuration for sim
            mat_config (dict): material configuration
            config_update (dict, optional): changes for geo_config,
                sim_config, mat_config. Defaults to {}.
            T_tp (float, optional): Triple point temperature. Defaults
                to 505 (tin).
            h_min (float, optional): Minimum allowable crystal height.
                Defaults to 0.001.
            h_max (float, optional): Maximum allowable crystal height.
                Defaults to 0.02.
            max_iterations (int, optional): Maximum number of
               iterations. Defaults to 10.
            dT_max (float, optional): Convergence criterion: maximum
                allowable difference in triple point temperature.
                Defaults to 0.01.
            sim_name (str, optional): simulation name.
                Defaults to "opencgs-simulation".
            base_dir (str, optional): path of base directory. Defaults
                to "./simdata".
            with_date (bool, optional): Include date in simulation
                directory name. Defaults to True.
            metadata (str, optional): Metadata to be saved, e.g git-hash
                of parent repository. Defaults to "".
        """
        super().__init__(
            geo,
            geo_config,
            sim,
            sim_config,
            mat_config,
            config_update,
            sim_name,
            "di",
            base_dir,
            with_date,
            metadata,
        )
        with open(self.input_dir + "/di_params.yml", "w") as f:
            yaml.dump(
                {
                    "T_tp": T_tp,
                    "h_min": h_min,
                    "h_max": h_max,
                    "max_iterations": max_iterations,
                    "dT_max": dT_max,
                },
                f,
            )
        self.T_tp = T_tp
        self.h_min = h_min
        self.h_max = h_max
        self.max_iterations = max_iterations
        self.dT_max = dT_max

    def execute(self):
        """Run iterative diameter computation process, consisting of
        multiple steady-state simulations."""
        # initial simulations
        geo_config = deepcopy(self.geo_config)
        geo_config["crystal"]["h"] = self.h_min
        sim_h_min = SteadyStateSim(
            self.geo,
            geo_config,
            self.sim,
            self.sim_config,
            self.mat_config,
            sim_name=f"h={self.h_min}",
            base_dir=self.sim_dir,
            with_date=False,
        )
        sim_h_min.execute()
        T_hmin = sim_h_min.T_tp
        geo_config = deepcopy(self.geo_config)
        geo_config["crystal"]["h"] = self.h_max
        sim_h_max = SteadyStateSim(
            self.geo,
            geo_config,
            self.sim,
            self.sim_config,
            self.mat_config,
            sim_name=f"h={self.h_max}",
            base_dir=self.sim_dir,
            with_date=False,
        )
        sim_h_max.execute()
        T_hmax = sim_h_max.T_tp
        # evaluate
        print(f"r-min: {self.h_min} m - T = {T_hmin:5f} K")
        print(f"r-max: {self.h_max} m - T = {T_hmax:5f} K")
        Ttp_h = {}  # triple point temperature : height
        Ttp_h.update({T_hmin: self.h_min})
        Ttp_h.update({T_hmax: self.h_max})
        self.export_Ttp_h(Ttp_h)
        if not T_hmin <= self.T_tp <= T_hmax:
            print("Diameter fitting impossible with current setup.")
            print("Interpolated height would be", self.compute_new_h(Ttp_h), "m.")
            with open(f"{self.res_dir}/iteration-summary.yml", "w") as f:
                yaml.dump(self.T_tp, f)
        # iteration
        converged = False
        for i in range(self.max_iterations):
            print("Diameter iteration", i + 1)
            h_new = self.compute_new_h(Ttp_h)
            print("new height:", h_new)
            if not (self.h_min < h_new < self.h_max):
                print("ERROR: Interpolation not possible.")
                break
            geo_config = deepcopy(self.geo_config)
            geo_config["crystal"]["h"] = h_new
            sim = SteadyStateSim(
                self.geo,
                geo_config,
                self.sim,
                self.sim_config,
                self.mat_config,
                sim_name=f"h={h_new}",
                base_dir=self.sim_dir,
                with_date=False,
            )
            sim.execute()
            Ttp_new = sim.T_tp
            print("corresponding TP temperature:", Ttp_new)
            Ttp_h.update({Ttp_new: h_new})
            self.export_Ttp_h(Ttp_h)
            if np.abs(Ttp_new - self.T_tp) <= self.dT_max:
                print("Iteration finished.")
                print("Crystal height =", h_new, "m.")
                print("TP Temperature =", Ttp_new, "K")
                converged = True
                break

        results = {
            "converged": converged,
            "iterations": i + 1,
            "dT at TP": Ttp_new - self.T_tp,
            "Ttp_h": Ttp_h,
        }
        with open(f"{self.res_dir}/iteration-summary.yml", "w") as f:
            yaml.dump(results, f)
        if converged:
            results.update({"height": h_new})
            shutil.copytree(sim.root_dir, f"{self.res_dir}/{sim.sim_name}")
        self.plot(Ttp_h)

    def plot(self, Ttp_h):
        """Plot convergence

        Args:
            Ttp_h (dict): Temperature @TP: crystal height
        """
        fig, ax = plt.subplots(1, 2, figsize=(5.75, 3))
        ax[0].plot(list(Ttp_h.keys()), "x-")
        ax[0].set_xlabel("simulation")
        ax[0].set_ylabel("temperature at triple point [K]")
        ax[0].grid()
        T_tps_sorted = {
            k: v for k, v in sorted(Ttp_h.items(), key=lambda item: item[1])
        }
        ax[1].plot(list(T_tps_sorted.values()), list(T_tps_sorted.keys()), "x-")
        ax[1].set_xlabel("crystal height [m]")
        ax[1].set_ylabel("temperature at triple point [K]")
        ax[1].grid()
        fig.tight_layout()
        fig.savefig(f"{self.res_dir}/diameter-iteration.png")
        plt.close(fig)

    def compute_new_h(self, Ttp_h):
        """Compute new crystal height for next iteration.

        Args:
            Ttp_h (dict): Temperature @TP: crystal height
                (from previous iterations)
        """
        T_tps = np.fromiter(Ttp_h.keys(), float)
        rs = np.fromiter(Ttp_h.values(), float)
        # ignore failed simulations (if possible)
        T_tps_new = np.array([T for T in T_tps if T > 0.0])
        if len(T_tps_new) >= 2:
            T_tps = T_tps_new
        # try to INTERpolate
        T1 = 0.0
        T2 = self.T_tp * 2
        for T in T_tps:
            if T1 < T < self.T_tp:
                T1 = T
            else:
                T2 = T
        try:
            h1 = Ttp_h[T1]
            h2 = Ttp_h[T2]
            h_new = h1 + (h2 - h1) / (T2 - T1) * (self.T_tp - T1)
            h_new = round(float(h_new), 8)
            if h_new == rs[-1]:
                print(
                    "WARNING: Non-linearity, tanking value from last iteration for interpolation."
                )
                if T_tps[-1] < self.T_tp:
                    T1 = T_tps[-1]
                    h1 = Ttp_h[T1]
                else:
                    T2 = T_tps[-1]
                    h2 = Ttp_h[T2]
        except KeyError:
            print("Warning: Could not interpolate!")
        # EXTRApolate if necessary
        if T1 == 0.0:
            T_diff = sorted(T_tps - self.T_tp)
            T1 = T_diff[0]
            T2 = T_diff[1]
        if T2 == self.T_tp * 2:
            T_diff = sorted(T_tps - self.T_tp)
            T1 = T_diff[-2]
            T2 = T_diff[-1]
        h_new = h1 + (h2 - h1) / (T2 - T1) * (self.T_tp - T1)
        h_new = round(float(h_new), 8)
        print("selected points for interpolation")
        print(Ttp_h)
        print("T1 =", T1)
        print("T2 =", T2)
        return h_new

    def export_Ttp_h(self, Ttp_h):
        """Write dictionary with triple point temperature and height to
        file."""
        with open(f"{self.res_dir}/Ttp_h.yml", "w") as f:
            yaml.dump(Ttp_h, f, sort_keys=False)
