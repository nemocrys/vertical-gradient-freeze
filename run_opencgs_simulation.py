# This function executes the simulation trough opencgs.
# Alternatively, you can run it directly with pyelmer.execute,
# see simulation_setup_pyelmer.py or simulation_setup_opencgs.py.

import opencgs.control as ctrl
from opencgs.sim import ParameterStudy, SteadyStateSim
from geometry import geometry
from simulation_setup_opencgs import simulation_opencgs


if __name__ == "__main__":
    try:
        git_metadata = ctrl.get_git_metadata()
    except:
        git_metadata = "not available"

    config_geo = ctrl.load_config("./config_geo.yml")
    config_sim = ctrl.load_config("./config_sim.yml")
    config_mat = ctrl.load_config("./config_mat.yml")
    config_opencgs = ctrl.load_config("./config_opencgs.yml")
    config_opencgs.update({"metadata": git_metadata})

    # This is used to run steady state simulations / parameter studies
    if "study_params" in config_opencgs:
        sim = ParameterStudy(
            SteadyStateSim,
            geometry,
            config_geo,
            simulation_opencgs,
            config_sim,
            config_mat,
            base_dir="simdata",
            **config_opencgs
        )
    else:
        sim = SteadyStateSim(
            geometry,
            config_geo,
            simulation_opencgs,
            config_sim,
            config_mat,
            base_dir="simdata",
            **config_opencgs
        )

    sim.execute()
