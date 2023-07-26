# implementation using pyelmer only
import os
from pyelmer import elmerkw as elmer
from pyelmer.execute import run_elmer_solver, run_elmer_grid
from pyelmer.post import scan_logfile
import yaml

from geometry import geometry


def simulation_pyelmer(
    model, config, sim_dir="./simdata", config_mat={}, elmer_config_file="config_elmer.yml"
):
    # setup simulation, solvers
    sim = elmer.load_simulation("axisymmetric_steady", elmer_config_file)
    solver_heat = elmer.load_solver("HeatSolver", sim, elmer_config_file)
    solver_phase_change = elmer.load_solver("SteadyPhaseChange", sim, elmer_config_file)
    solver_mesh = elmer.load_solver("MeshUpdate", sim, elmer_config_file)
    solver_out = elmer.load_solver("ResultOutputSolver", sim, elmer_config_file)
    eqn_main = elmer.Equation(sim, "eqn_main", [solver_heat, solver_mesh])
    eqn_phase_change = elmer.Equation(sim, "eqn_phase_change", [solver_phase_change])

    # add crystal
    crystal = elmer.Body(sim, "crystal", [model["crystal"].ph_id])
    material_name = model["crystal"].params.material
    mat = elmer.Material(sim, material_name, config_mat[material_name])
    melting_point = mat.data["Melting Point"]
    ic = elmer.InitialCondition(
        sim, "T_crystal", {"Temperature": model["crystal"].params.T_init}
    )
    crystal.equation = eqn_main
    crystal.material = mat
    crystal.initial_condition = ic

    # add melt
    melt = elmer.Body(sim, "melt", [model["melt"].ph_id])
    material_name = model["melt"].params.material
    mat = elmer.Material(sim, material_name, config_mat[material_name])
    ic = elmer.InitialCondition(
        sim, "T_melt", {"Temperature": model["melt"].params.T_init}
    )
    melt.equation = eqn_main
    melt.material = mat
    melt.initial_condition = ic

    # add heaters
    for shape in [
        "heater_side_bot",
        "heater_side_mid",
        "heater_side_top",
        "heater_top",
    ]:
        force = elmer.BodyForce(
            sim,
            "force_" + shape,
            {
                "Heat Source": 1,
                "Integral Heat Source": config["heating"][shape],
                "Smart Heater Control": "Logical True",
            },
        )
        bdy = elmer.Body(sim, shape, [model[shape].ph_id])
        material_name = model[shape].params.material
        mat = elmer.Material(sim, material_name, config_mat[material_name])
        ic = elmer.InitialCondition(
            sim, "T_" + shape, {"Temperature": model[shape].params.T_init}
        )
        bdy.equation = eqn_main
        bdy.material = mat
        bdy.body_force = force
        bdy.initial_condition = ic

    # add other bodies
    for shape in [
        "base_plate",
        "insulation_bot",
        "crucible",
        "enclosure",
        "insulation",
    ]:
        bdy = elmer.Body(sim, shape, [model[shape].ph_id])
        material_name = model[shape].params.material
        mat = elmer.Material(sim, material_name, config_mat[material_name])
        bdy.equation = eqn_main
        bdy.material = mat

    # setup phase change
    t0_phase_change = elmer.InitialCondition(
        sim, "t0_phase_change", {"Temperature": melting_point}
    )
    melt_crystal_if = elmer.Body(
        sim, "melt_crystal_if", [model["if_melt_crystal"].ph_id]
    )
    melt_crystal_if.equation = eqn_phase_change
    melt_crystal_if.material = crystal.material
    melt_crystal_if.initial_condition = t0_phase_change

    if_melt_crystal = elmer.Boundary(
        sim,
        "if_melt_crystal",
        [model["if_melt_crystal"].ph_id],
    )
    if_melt_crystal.smart_heater = True
    if_melt_crystal.smart_heater_T = melting_point
    if_melt_crystal.phase_change_steady = True
    if_melt_crystal.phase_change_body = melt_crystal_if
    if_melt_crystal.normal_target_body = crystal
    if_melt_crystal.material = crystal.material

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
        bnd = elmer.Boundary(sim, bnd, [model[bnd].ph_id])
        bnd.radiation = True
        bnd.mesh_update = [0, 0]

    # add outside boundaries
    for bnd in [
        "bnd_baseplate_outside",
        "bnd_insulation_outside",
    ]:
        bnd = elmer.Boundary(sim, bnd, [model[bnd].ph_id])
        bnd.fixed_temperature = config["T_ambient"]
        bnd.mesh_update = [0, 0]

    # add interfaces
    bnd = elmer.Boundary(sim, "symmetry_axis", [model["symmetry_axis"].ph_id])
    bnd.mesh_update = [0, None]

    for bnd in [
        "if_crucible_melt",
        "if_crucible_crystal",
        "if_crucible_insbot",
        "if_insbot_enclosure",
        "if_base_enclosure",
        "if_base_insulation",
    ]:
        bnd = elmer.Boundary(sim, bnd, [model[bnd].ph_id])
        bnd.mesh_update = [0, 0]

    # write Elmer simulation input file (sif-file)
    sim.write_sif(sim_dir)


if __name__ == "__main__":
    sim_dir = "./simdata/pyelmer_simulation"
    if os.path.exists(sim_dir):
        raise ValueError("Please remove the old simulation directory.")

    with open("config_geo.yml") as f:
        config_geo = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir)

    with open("config_sim.yml") as f:
        config_sim = yaml.safe_load(f)
    with open("config_mat.yml") as f:
        config_mat = yaml.safe_load(f)

    simulation_pyelmer(model, config_sim, sim_dir, config_mat)

    run_elmer_grid(sim_dir, "case.msh")
    run_elmer_solver(sim_dir)
    err, warn, stats = scan_logfile(sim_dir)
    print("Errors:", err)
    print("Warnings:", warn)
    print("Statistics:", stats)
