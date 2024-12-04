# This function executes the simulation trough opencgs.
# Alternatively, you can run it directly with pyelmer.execute,
# see simulation_setup_pyelmer.py or simulation_setup_opencgs.py.

import opencgs.control as ctrl
from opencgs.sim import ParameterStudy
from geometry import geometry
from simulation_setup_opencgs_fixed_temperature import simulation_opencgs

from diameter_iteraion_vgf import DiameterIterationVGF


if __name__ == "__main__":
    try:
        git_metadata = ctrl.get_git_metadata()
    except:
        git_metadata = "not available"

    config_geo = ctrl.load_config("./config_geo.yml")
    config_sim = ctrl.load_config("config_sim_fixed-temperature.yml")
    config_mat = ctrl.load_config("./config_mat.yml")
    config_opencgs = ctrl.load_config("./config_opencgs_fixed-temperature.yml")
    config_opencgs.update({"metadata": git_metadata})

    # This is used to run steady state simulations / parameter studies
    # TODO this gives an error due to a failing type check. Should be easily fixed.
    if "study_params" in config_opencgs:
        sim = ParameterStudy(
            DiameterIterationVGF,
            geometry,
            config_geo,
            simulation_opencgs,
            config_sim,
            config_mat,
            base_dir="simdata",
            **config_opencgs
        )
    else:
        sim = DiameterIterationVGF(
            geometry,
            config_geo,
            simulation_opencgs,
            config_sim,
            config_mat,
            base_dir="simdata",
            **config_opencgs
        )

    sim.execute()
