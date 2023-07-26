import os
import numpy as np
from objectgmsh import Model, Shape, MeshControlLinear, MeshControlExponential, cut
import gmsh
import yaml

from pyelmer import elmerkw as elmer
from opencgs.setup import ElmerSetupCz


occ = gmsh.model.occ


def geometry(config, sim_dir="./simdata", visualize=False):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    model = Model()

    # base plate is a cylinder in 3D = a rectangle in 2D axisymmetric
    base_plate = Shape(
        model,
        2,
        "base_plate",
        [
            occ.add_rectangle(
                0,
                -config["base_plate"]["t"],
                0,
                config["base_plate"]["r"],
                config["base_plate"]["t"],
            )
        ],
    )
    base_plate.mesh_size = config["base_plate"]["mesh_size"]
    base_plate.params.material = config["base_plate"][
        "material"
    ]  # we use "params" to save various values we need later

    # bottom insulation is a rectangle
    insulation_bot = Shape(
        model,
        2,
        "insulation_bot",
        [
            occ.add_rectangle(
                0,
                config["insulation_bot"]["y0"],
                0,
                config["insulation_bot"]["r"],
                config["insulation_bot"]["h"],
            )
        ],
    )
    insulation_bot.mesh_size = config["insulation_bot"]["mesh_size"]
    insulation_bot.params.material = config["insulation_bot"]["material"]
    insulation_bot_y_top = (
        config["insulation_bot"]["y0"] + config["insulation_bot"]["h"]
    )

    # crucible is a rectangle with a hole (used for melt later)
    crucible = Shape(
        model,
        2,
        "crucible",
        [
            occ.add_rectangle(
                0,
                insulation_bot_y_top,
                0,
                config["crucible"]["r_out"],
                config["crucible"]["h"],
            )
        ],
    )
    crucible.mesh_size = config["crucible"]["mesh_size"]
    crucible.params.material = config["crucible"]["material"]
    crucible.params.T_init = config["crucible"]["T_init"]
    crucible_y_in = insulation_bot_y_top + config["crucible"]["t_bot"]
    crucible_h_cone = config["crucible"]["r_in"] / np.tan(
        np.deg2rad(config["crucible"]["angle"] / 2)
    )
    points = [
        occ.add_point(0, crucible_y_in, 0),
        occ.add_point(
            config["crucible"]["r_in"],
            crucible_y_in + crucible_h_cone,
            0,
        ),
        occ.add_point(
            config["crucible"]["r_in"],
            insulation_bot_y_top + config["crucible"]["h"],
            0,
        ),
        occ.add_point(
            0,
            insulation_bot_y_top + config["crucible"]["h"],
            0,
        ),
    ]
    lines = [occ.add_line(points[i - 1], points[i]) for i in range(len(points))]
    loop = occ.add_curve_loop(lines)
    melt_volume = occ.add_surface_filling(loop)
    occ.cut(crucible.dimtags, [(2, melt_volume)], removeTool=False)

    # crystal
    crystal = Shape(
        model,
        2,
        "crystal",
    )
    crystal.mesh_size = config["crystal"]["mesh_size"]
    crystal.params.material = config["crystal"]["material"]
    crystal.params.T_init = config["crystal"]["T_init"]
    crystal_y_top = crucible_y_in + config["crystal"]["h"]
    points = [
        occ.add_point(0, crucible_y_in, 0),
    ]
    if config["crystal"]["h"] < crucible_h_cone:
        points.append(
            occ.add_point(
                config["crystal"]["h"]
                * np.tan(np.deg2rad(config["crucible"]["angle"] / 2)),
                crucible_y_in + config["crystal"]["h"],
                0,
            )
        )
    else:
        points.extend(
            [
                occ.add_point(
                    config["crucible"]["r_in"], crucible_y_in + crucible_h_cone, 0
                ),
                occ.add_point(config["crucible"]["r_in"], crystal_y_top, 0),
            ]
        )
    points.append(occ.add_point(0, crystal_y_top, 0))
    lines = [occ.add_line(points[i - 1], points[i]) for i in range(len(points))]
    loop = occ.add_curve_loop(lines)
    crystal.geo_ids = [occ.add_surface_filling(loop)]

    # melt
    melt = Shape(
        model, 2, "melt", [melt_volume]
    )  # we re-use the part that we used to cut the hole into the crucible
    melt.mesh_size = config["melt"]["mesh_size"]
    melt.params.material = config["melt"]["material"]
    melt.params.T_init = config["melt"]["T_init"]
    occ.cut(melt.dimtags, crystal.dimtags, removeTool=False)
    cutbox = occ.add_rectangle(
        0,
        crystal_y_top + config["melt"]["h"],
        0,
        config["crucible"]["r_in"],
        config["crucible"]["h"],
    )
    occ.cut(melt.dimtags, [(2, cutbox)])

    # enclosure
    enclosure = Shape(
        model,
        2,
        "enclosure",
        [
            occ.add_rectangle(
                0,
                0,
                0,
                config["enclosure"]["r_in"] + config["enclosure"]["t"],
                config["enclosure"]["h_in"] + config["enclosure"]["t"] * 2,
            )
        ],
    )
    enclosure.mesh_size = config["enclosure"]["mesh_size"]
    enclosure.params.material = config["enclosure"]["material"]
    hole = occ.add_rectangle(
        0,
        config["enclosure"]["t"],
        0,
        config["enclosure"]["r_in"],
        config["enclosure"]["h_in"],
    )
    occ.cut(enclosure.dimtags, [(2, hole)])

    # heaters
    heaters = []
    for htr in ["heater_side_bot", "heater_side_mid", "heater_side_top", "heater_top"]:
        heater = Shape(
            model,
            2,
            htr,
            [
                occ.add_rectangle(
                    config[htr]["r_in"],
                    config[htr]["y0"],
                    0,
                    config[htr]["dx"],
                    config[htr]["dy"],
                )
            ],
        )
        heater.mesh_size = config[htr]["mesh_size"]
        heater.params.material = config[htr]["material"]
        heater.params.T_init = config[htr]["T_init"]
        heaters.append(heater)

    # insulation
    insulation = Shape(
        model,
        2,
        "insulation",
        [
            occ.add_rectangle(
                0,
                0,
                0,
                config["insulation"]["r_in"] + config["insulation"]["t"],
                config["insulation"]["h_in"] + config["insulation"]["t"],
            )
        ],
    )
    insulation.mesh_size = config["insulation"]["mesh_size"]
    insulation.params.material = config["insulation"]["material"]
    atmosphere = occ.add_rectangle(
        0, 0, 0, config["insulation"]["r_in"], config["insulation"]["h_in"]
    )  # will be used bellow
    occ.cut(insulation.dimtags, [(2, atmosphere)], removeTool=False)

    # atmosphere, will be used to determine surfaces
    # it is just a helper shape and will be removed later
    shapes = model.get_shapes(2)
    atmosphere = Shape(model, 2, "atmosphere", [atmosphere])
    for shape in shapes:
        atmosphere.geo_ids = cut(atmosphere.dimtags, shape.dimtags, remove_tool=False)

    # set interfaces between shapes, this removes duplicate lines and ensures a consistent mesh
    occ.fragment(
        base_plate.dimtags
        + insulation_bot.dimtags
        + crucible.dimtags
        + crystal.dimtags
        + melt.dimtags
        + enclosure.dimtags
        + atmosphere.dimtags,
        [],
    )

    model.synchronize()

    # extract phase interface
    if_melt_crystal = Shape(model, 1, "if_melt_crystal", melt.get_interface(crystal))

    # extract boundaries for surface-to-surface radiation
    bnd_baseplate = Shape(
        model, 1, "bnd_baseplate", base_plate.get_interface(atmosphere)
    )
    bnd_insulation_bot = Shape(
        model, 1, "bnd_insulation_bot", insulation_bot.get_interface(atmosphere)
    )
    bnd_crucible = Shape(model, 1, "bnd_crucible", crucible.get_interface(atmosphere))
    bnd_melt = Shape(model, 1, "bnd_melt", melt.get_interface(atmosphere))
    bnd_enclosure = Shape(
        model, 1, "bnd_enclosure", enclosure.get_interface(atmosphere)
    )
    for heater in heaters:
        bnd = Shape(model, 1, "bnd_" + heater.name, heater.get_interface(atmosphere))
    bnd_insulation = Shape(
        model, 1, "bnd_insulation", insulation.get_interface(atmosphere)
    )

    bnd_baseplate_outside = Shape(
        model,
        1,
        "bnd_baseplate_outside",
        [base_plate.bottom_boundary, base_plate.right_boundary],
    )
    bnd_insulation_outside = Shape(
        model,
        1,
        "bnd_insulation_outside",
        [insulation.right_boundary, insulation.top_boundary],
    )

    # extract interfaces (not all are required but gives much better visualization in ParaView)
    if_crucible_melt = Shape(
        model,
        1,
        "if_crucible_melt",
        crucible.get_interface(melt),
    )
    if_crucible_crystal = Shape(
        model,
        1,
        "if_crucible_crystal",
        crucible.get_interface(crystal),
    )
    if_crucible_insbot = Shape(
        model,
        1,
        "if_crucible_insbot",
        crucible.get_interface(insulation_bot),
    )
    if_insbot_enclosure = Shape(
        model,
        1,
        "if_insbot_enclosure",
        enclosure.get_interface(insulation_bot),
    )
    if_base_enclosure = Shape(
        model,
        1,
        "if_base_enclosure",
        enclosure.get_interface(base_plate),
    )
    if_base_insulation = Shape(
        model,
        1,
        "if_base_insulation",
        insulation.get_interface(base_plate),
    )

    # symmetry axis
    symmetry_axis = Shape(model, 1, "symmetry_axis", model.symmetry_axis)

    model.remove_shape(atmosphere)

    model.make_physical()

    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()
    model.generate_mesh(**config["mesh"])
    if visualize:
        model.show()
    print(model)
    model.write_msh(sim_dir + "/case.msh")

    model.close_gmsh()
    return model


def simulation_pyelmer(
    model, config, config_mat, sim_dir="./simdata", elmer_config_file="config_elmer.yml"
):
    # implementation using pyelmer only

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


def simulation_opencgs(model, config, config_mat, sim_dir="./simdata"):
    # alternative implementation using opencgs czochralski simulation
    sim = ElmerSetupCz(
        sim_dir=sim_dir,
        heat_control=True,
        heat_convection=False,
        phase_change=True,
        heating_resistance=True,
        smart_heater={
            "T": 505,  # melting point Sn
            "control-point": False,  # if False: use T at melt-crystal interface (triple point)
            # "x": 0.035,
            # "y": 0.005,
            # "z": 0.0,
        },
        solver_update={
            "global": {"Steady State Max Iterations": 10},
            "all-solvers": {
                "Linear System Iterative Method": "Idrs",
                "Linear System Residual Output": 10,
            },
        },
        probes={"seed": [0.0, 0.0]},
        materials_dict=config_mat,
    )
    sim.sim.solvers["ResultOutputSolver"].data.update({"Vtu Part collection": True})

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
        heater.data.update({"name": shape})

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
    sim_dir = "./simdata"

    with open("config_geo.yml") as f:
        config_geo = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir)

    with open("config_sim.yml") as f:
        config_sim = yaml.safe_load(f)
    with open("config_mat.yml") as f:
        config_mat = yaml.safe_load(f)

    # simulation_pyelmer(model, config_sim, config_mat, sim_dir)
    simulation_opencgs(model, config_sim, config_mat, sim_dir)

    from pyelmer.execute import run_elmer_solver, run_elmer_grid
    from pyelmer.post import scan_logfile

    run_elmer_grid(sim_dir, "case.msh")
    run_elmer_solver(sim_dir)
    err, warn, stats = scan_logfile(sim_dir)
    print("Errors:", err)
    print("Warnings:", warn)
    print("Statistics:", stats)
