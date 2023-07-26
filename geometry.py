import os
import numpy as np
from objectgmsh import Model, Shape, MeshControlLinear, MeshControlExponential, cut
import gmsh
import yaml


occ = gmsh.model.occ


def geometry(config, sim_dir="./simdata", name="vgf", visualize=False):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    model = Model(name)

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
        crystal_y_top + config["melt"]["h"] - config["crystal"]["h"],
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

    # mesh settings
    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()

    # add linear mesh control to ensure smooth transition in mesh sizes
    max_meshsize = 0.1
    for shape in model.get_shapes(2):
        MeshControlLinear(model, shape, shape.mesh_size, max_meshsize)

    # add refinement at crystal melt interface
    MeshControlExponential(model, if_melt_crystal, melt.mesh_size / 5, exp=1.4)

    model.generate_mesh(**config["mesh"])
    if visualize:
        model.show()
    print(model)
    model.write_msh(sim_dir + "/case.msh")

    model.close_gmsh()
    return model


if __name__ == "__main__":
    sim_dir = "./simdata"

    with open("config_geo.yml") as f:
        config_geo = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir, visualize=True)
