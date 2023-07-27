# alternative implementation using opencgs czochralski simulation
from opencgs.setup import ElmerSetupCz
import os
from pyelmer.execute import run_elmer_solver, run_elmer_grid
from pyelmer.post import scan_logfile
import yaml

from geometry import geometry


def simulation_opencgs(model, config, sim_dir="./simdata", config_mat={}):
    sim = ElmerSetupCz(
        sim_dir=sim_dir,
        heat_control=True,
        heat_convection=False,
        phase_change=True,
        heating_resistance=True,
        smart_heater={
            "T": 1511,  # melting point GaAs
            "control-point": False,  # if False: use T at melt-crystal interface (triple point)
            # "x": 0.035,
            # "y": 0.005,
            # "z": 0.0,
        },
        solver_update={  # some changes to the default opencgs solver settings to make simulation faster
            "global": {"Steady State Max Iterations": 10},
            "all-solvers": {
                "Linear System Iterative Method": "Idrs",
                "Linear System Residual Output": 10,
            },
        },
        probes={  # coordinates where temperature is evaluated
            "seed": [0.0, 0.152],
            "heater_side_bot": [0.109, 0.122],
            "heater_side_mid": [0.109, 0.276],
            "heater_side_top": [0.109, 0.43],
            "heater_top": [0.0, 0.509],
            },
        materials_dict=config_mat,
    )
    sim.sim.solvers["ResultOutputSolver"].data.update({"Vtu Part collection": True})  # gives nicer output

    # add crystal
    crystal = sim.add_crystal(model["crystal"])
    crystal.data.update({"name": "crystal"})

    # add heaters
    for shape in [
        "heater_side_bot",
        "heater_side_mid",
        "heater_side_top",
        "heater_top",
    ]:
        heater = sim.add_resistance_heater(model[shape], config["heating"][shape])
        heater.data.update({"name": shape})  # gives nicer output

    # add other bodies
    for body in [
        "base_plate",
        "insulation_bot",
        "crucible",
        "melt",
        "enclosure",
        "insulation",
    ]:
        bdy = sim.add_body(model[body])
        bdy.data.update({"name": body})

    # phase interface
    sim.add_phase_interface(model["if_melt_crystal"])

    # add boundaries with surface-to-surface radiation
    for bnd in [
        "bnd_baseplate",
        "bnd_insulation_bot",
        "bnd_crucible",
        "bnd_melt",
        "bnd_enclosure",
        "bnd_heater_side_bot",
        "bnd_heater_side_mid",
        "bnd_heater_side_top",
        "bnd_heater_top",
        "bnd_insulation",
    ]:
        sim.add_radiation_boundary(model[bnd])

    # add outside boundaries
    for bnd in [
        "bnd_baseplate_outside",
        "bnd_insulation_outside",
    ]:
        sim.add_temperature_boundary(model[bnd], config["T_ambient"])

    # interfaces
    sim.add_interface(
        model["symmetry_axis"], movement=[0, None]
    )  # allow deformation in axial direction
    for interface in [
        "if_crucible_melt",
        "if_crucible_crystal",
        "if_crucible_insbot",
        "if_insbot_enclosure",
        "if_base_enclosure",
        "if_base_insulation",
    ]:
        sim.add_interface(model[interface])
    # export
    sim.export()


if __name__ == "__main__":
    # Option 1: run directly with pyelmer.execute (see below)
    # Option 2: run with opencgs (see run_opencgs.py)
    sim_dir = "./simdata/opencgs_simulation"
    if os.path.exists(sim_dir):
        raise ValueError("Please remove the old simulation directory.")

    with open("config_geo.yml") as f:
        config_geo = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir)

    with open("config_sim.yml") as f:
        config_sim = yaml.safe_load(f)
    with open("config_mat.yml") as f:
        config_mat = yaml.safe_load(f)

    simulation_opencgs(model, config_sim, sim_dir, config_mat)

    run_elmer_grid(sim_dir, "case.msh")
    run_elmer_solver(sim_dir)
    err, warn, stats = scan_logfile(sim_dir)
    print("Errors:", err)
    print("Warnings:", warn)
    print("Statistics:", stats)
