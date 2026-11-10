from collections import deque
from math import ceil
from typing import Any
from vpython import color, curve, cylinder, rate, scene, sphere, vector

from dirabo.chem.i_keybd import bind_keyboard_controls
from dirabo.chem.math import Vec3, add, length, scale

from dirabo.chem.particles import (
    SIMULATION_DT,
    SIMULATION_STEPS_PER_FRAME,
    NUCLEUS_REPULSION_FACTOR,
    TRACE_MARKER_INTERVAL_SECONDS,
    TRACE_POINT_LIMIT,
    Electron,
    Molecule,
    MoleculeCluster,
    NuclearTrace,
    _force_color,
    _forces_by_particle,
    integrate_electrons,
    _particle_states,
    recenter_system,
    force_report,
)


def _record_nuclear_trace(
    trace: NuclearTrace,
    position: Vec3,
    sample_index: int,
    marker_interval: int,
    marker_capacity: int,
) -> None:
    trace.positions.append(position)
    trace.visual.append(pos=vector(*position))

    if sample_index % TRACE_POINT_LIMIT == 0:
        trace.visual.clear()
        for saved_position in trace.positions:
            trace.visual.append(pos=vector(*saved_position))

    if sample_index % marker_interval != 0:
        return

    marker_size = 0.035
    if len(trace.marker_visuals) < marker_capacity:
        x_marker = cylinder(
            pos=vector(position[0] - marker_size / 2, position[1], position[2]),
            axis=vector(marker_size, 0, 0),
            radius=0.004,
            color=color.cyan,
        )
        y_marker = cylinder(
            pos=vector(position[0], position[1] - marker_size / 2, position[2]),
            axis=vector(0, marker_size, 0),
            radius=0.004,
            color=color.cyan,
        )
        x_marker.visible = trace.visual.visible
        y_marker.visible = trace.visual.visible
        trace.marker_visuals.append((x_marker, y_marker))
        return

    x_marker, y_marker = trace.marker_visuals[trace.next_marker_slot]
    x_marker.pos = vector(position[0] - marker_size / 2, position[1], position[2])
    x_marker.axis = vector(marker_size, 0, 0)
    y_marker.pos = vector(position[0], position[1] - marker_size / 2, position[2])
    y_marker.axis = vector(0, marker_size, 0)
    trace.next_marker_slot = (trace.next_marker_slot + 1) % marker_capacity


def show_molecule(molecule: Molecule) -> None:
    molecules = MoleculeCluster(molecule)
    show_molecules(molecules)


def show_molecules(molecules: MoleculeCluster) -> None:
    atoms = molecules.atoms
    electrons = molecules.electrons
    bonds = molecules.bonds
    lone_pairs = molecules.lone_pairs
    if TRACE_POINT_LIMIT < 2:
        raise ValueError("TRACE_POINT_LIMIT muss mindestens 2 sein")
    if TRACE_MARKER_INTERVAL_SECONDS <= 0:
        raise ValueError("Das Trace-Marker-Intervall muss positiv sein")
    simulation_frame_dt = SIMULATION_DT * SIMULATION_STEPS_PER_FRAME
    if simulation_frame_dt <= 0:
        raise ValueError("Die Simulationszeit pro Bild muss positiv sein")
    points_per_second = 1 / simulation_frame_dt
    marker_interval = max(
        1,
        round(points_per_second * TRACE_MARKER_INTERVAL_SECONDS),
    )
    marker_capacity = ceil(TRACE_POINT_LIMIT / marker_interval) + 1
    atom_ids = [atom.atom_id for atom in atoms]
    if len(atom_ids) != len(set(atom_ids)):
        raise ValueError("Atom-IDs müssen über alle Moleküle hinweg eindeutig sein")

    scene.title = ", ".join(molecule.name for molecule in molecules.molecules)
    scene.width = 1000
    scene.height = 700
    scene.background = color.black
    scene.range = max(1.6, 1.6 * len(molecules.molecules))

    atoms_by_id = {atom.atom_id: atom for atom in atoms}
    for electron in electrons:
        electron.initialize(atoms_by_id[electron.atom_id])
    recenter_system(molecules, remove_com_velocity=False)

    traces = {
        atom.atom_id: NuclearTrace(
            positions=deque([atom.position], maxlen=TRACE_POINT_LIMIT),
            visual=curve(
                pos=[vector(*atom.position)],
                color=color.cyan,
                radius=0.004,
                retain=TRACE_POINT_LIMIT,
            ),
        )
        for atom in atoms
    }

    for atom in atoms:
        atom_color = color.red if atom.element == "O" else color.white
        atom.visual = sphere(
            pos=vector(*atom.position),
            radius=atom.display_radius,
            color=atom_color,
        )
        atom.normal_display_radius = atom.display_radius
        atom.force_visual = cylinder(
            pos=vector(*atom.position),
            axis=vector(0, 0, 0),
            radius=0.008,
            color=color.red,
        )

    for bond in bonds:
        bond.visual = cylinder(
            pos=vector(0, 0, 0),
            axis=vector(0, 0, 0),
            radius=0.035,
            color=color.white,
        )

    electrons_by_id: dict[str, Electron] = {}
    for electron in electrons:
        electrons_by_id[electron.electron_id] = electron
        electron.visual = sphere(
            pos=vector(*electron.position),
            radius=electron.display_radius,
            color=color.yellow,
        )
        electron.force_visual = cylinder(
            pos=vector(*electron.position),
            axis=vector(0, 0, 0),
            radius=0.008,
            color=color.red,
        )

    for target in electrons_by_id.values():
        for source in electrons_by_id.values():
            if source.electron_id == target.electron_id:
                continue
            target.individual_force_visuals[source.electron_id] = cylinder(
                pos=vector(*target.position),
                axis=vector(0, 0, 0),
                radius=0.004,
                color=color.blue,
            )

    for pair in lone_pairs:
        pair.visual = cylinder(
            pos=vector(0, 0, 0),
            axis=vector(0, 0, 0),
            radius=0.035,
            color=vector(0.82, 0.75, 0.62),
        )

    nuclear_force_visuals: dict[tuple[str, ...], tuple[Any, Any]] = {}
    selected_force_visuals: dict[
        tuple[str, ...],
        tuple[Any, str, str, float, Any, float, Vec3],
    ] = {}
    controls = bind_keyboard_controls(
        molecules,
        electrons_by_id,
        traces,
        nuclear_force_visuals,
        selected_force_visuals,
    )

    trace_sample_index = 0
    while True:
        rate(60)
        if controls.paused:
            continue
        for _ in range(SIMULATION_STEPS_PER_FRAME):
            integrate_electrons(
                molecules,
                SIMULATION_DT,
                inter_molecular_bond_forces_enabled=(controls.inter_molecular_bond_forces_enabled),
            )

        trace_sample_index += 1
        for atom in atoms:
            atom.visual.pos = vector(*atom.position)
            _record_nuclear_trace(
                traces[atom.atom_id],
                atom.position,
                trace_sample_index,
                marker_interval,
                marker_capacity,
            )

        for electron in electrons_by_id.values():
            electron.visual.pos = vector(*electron.position)

        for bond in bonds:
            first, second = bond.electrons
            start = first.visual.pos
            end = second.visual.pos
            bond.visual.pos = start
            bond.visual.axis = end - start

        report = force_report(
            molecules,
            inter_molecular_bond_forces_enabled=(controls.inter_molecular_bond_forces_enabled),
        )
        forces_by_particle = _forces_by_particle(report)
        particle_positions = {particle.particle_id: particle.position for particle in _particle_states(molecules)}
        atoms_by_nucleus_id = {f"{atom.atom_id}:nucleus": atom for atom in atoms}
        electrons_by_id = {electron.electron_id: electron for electron in electrons}
        selected_force_visuals.clear()
        nuclear_force_components: dict[
            tuple[str, ...],
            tuple[Vec3, Vec3, str, Vec3],
        ] = {}
        for contribution in report.contributions:
            target_atom = atoms_by_nucleus_id.get(contribution.target_id)
            origin_kind = "direct"
            target_nucleus_id = contribution.target_id
            if target_atom is None:
                target_electron = electrons_by_id.get(contribution.target_id)
                if target_electron is None or target_electron.bond_id is None:
                    continue
                target_atom = next(atom for atom in atoms if atom.atom_id == target_electron.atom_id)
                target_nucleus_id = f"{target_atom.atom_id}:nucleus"
                origin_kind = "bond-reaction"

            source_position = particle_positions.get(contribution.source_id)
            if source_position is None:
                continue
            component_force = contribution.force
            if origin_kind == "direct" and contribution.interaction == "nucleus-nucleus":
                component_force = scale(component_force, NUCLEUS_REPULSION_FACTOR)
            key = (
                target_nucleus_id,
                contribution.source_id,
                contribution.interaction,
                origin_kind,
                contribution.behavior,
            )
            previous = nuclear_force_components.get(key)
            if previous is not None:
                component_force = add(previous[3], component_force)
            nuclear_force_components[key] = (
                source_position,
                target_atom.position,
                contribution.behavior,
                component_force,
            )

        nuclear_resultants: dict[str, Vec3] = {}
        for key, (_, _, _, component_force) in nuclear_force_components.items():
            nuclear_resultants[key[0]] = add(
                nuclear_resultants.get(key[0], (0.0, 0.0, 0.0)),
                component_force,
            )

        visible_nuclear_force_keys = set(nuclear_force_components)
        for key in tuple(nuclear_force_visuals):
            if key not in visible_nuclear_force_keys:
                source_line, force_vector = nuclear_force_visuals[key]
                source_line.visible = False
                force_vector.visible = False

        for key, (
            source_position,
            nucleus_position,
            behavior,
            component_force,
        ) in nuclear_force_components.items():
            visuals = nuclear_force_visuals.get(key)
            if visuals is None:
                source_line = cylinder(
                    pos=vector(*source_position),
                    axis=vector(*nucleus_position) - vector(*source_position),
                    radius=0.002,
                    color=color.gray(0.6),
                )
                force_vector = cylinder(
                    pos=vector(*nucleus_position),
                    axis=vector(0, 0, 0),
                    radius=0.006,
                    color=_force_color(behavior),
                )
                source_line.visible = controls.nuclear_force_traces_visible
                force_vector.visible = controls.nuclear_force_traces_visible
                nuclear_force_visuals[key] = (source_line, force_vector)
            else:
                source_line, force_vector = visuals
                source_line.visible = controls.nuclear_force_traces_visible
                force_vector.visible = controls.nuclear_force_traces_visible

            source_line.pos = vector(*source_position)
            source_line.axis = vector(*nucleus_position) - vector(*source_position)
            force_length = length(component_force)
            displayed_length = min(0.55, force_length * 4.0)
            if force_length > 0:
                displayed_force = scale(
                    component_force,
                    displayed_length / force_length,
                )
            else:
                displayed_force = (0.0, 0.0, 0.0)
            force_vector.pos = vector(*nucleus_position)
            force_vector.axis = vector(*displayed_force)
            original_color = _force_color(behavior)
            force_length = length(component_force)
            force_vector.color = (
                vector(0.56, 0.0, 1.0) if controls.selected_force_key == key and force_length > 0 else original_color
            )
            if controls.nuclear_force_traces_visible and force_length > 0:
                selected_force_visuals[key] = (
                    force_vector,
                    key[1],
                    key[2],
                    force_length,
                    original_color,
                    displayed_length,
                    component_force,
                )

        electron_force_components: dict[
            str,
            dict[str, tuple[Vec3, set[str], set[str]]],
        ] = {electron_id: {} for electron_id in electrons_by_id}
        for contribution in report.contributions:
            if contribution.target_id not in electrons_by_id:
                continue
            sources = electron_force_components[contribution.target_id]
            component, behaviors, interactions = sources.get(
                contribution.source_id,
                ((0.0, 0.0, 0.0), set(), set()),
            )
            sources[contribution.source_id] = (
                add(component, contribution.force),
                behaviors | {contribution.behavior},
                interactions | {contribution.interaction},
            )

        for electron in electrons_by_id.values():
            force = forces_by_particle.get(electron.electron_id, (0.0, 0.0, 0.0))
            force_length = length(force)
            displayed_length = min(0.85, force_length * 4.0)
            if controls.electron_forces_visible and force_length > 0:
                displayed_force = scale(force, displayed_length / force_length)
            else:
                displayed_force = (0.0, 0.0, 0.0)
            electron.force_visual.pos = electron.visual.pos
            electron.force_visual.axis = vector(*displayed_force)
            resultant_key = (
                electron.electron_id,
                "resultant",
                "resultant",
            )
            electron.force_visual.color = (
                vector(0.56, 0.0, 1.0)
                if controls.selected_force_key == resultant_key and force_length > 0
                else color.yellow
            )
            sources = electron_force_components[electron.electron_id]
            for source_id, (individual_force, behaviors, interactions) in sources.items():
                force_length = length(individual_force)
                displayed_length = min(0.55, force_length * 4.0)
                if force_length > 0:
                    displayed_force = scale(individual_force, displayed_length / force_length)
                else:
                    displayed_force = (0.0, 0.0, 0.0)

                force_visual = electron.individual_force_visuals.get(source_id)
                if force_visual is None:
                    force_visual = cylinder(
                        pos=electron.visual.pos,
                        axis=vector(0, 0, 0),
                        radius=0.006,
                        color=color.gray(0.6),
                    )
                    force_visual.visible = controls.electron_forces_visible
                    electron.individual_force_visuals[source_id] = force_visual
                force_visual.visible = controls.electron_forces_visible
                force_visual.pos = electron.visual.pos
                force_visual.axis = vector(*displayed_force)
                behavior = next(iter(behaviors)) if len(behaviors) == 1 else "balanced"
                original_color = _force_color(behavior)
                force_key = (
                    electron.electron_id,
                    source_id,
                    "particle-force",
                )
                force_visual.color = (
                    vector(0.56, 0.0, 1.0)
                    if controls.selected_force_key == force_key and force_length > 0
                    else original_color
                )
                if controls.electron_forces_visible and force_length > 0:
                    interaction = " + ".join(sorted(interactions))
                    selected_force_visuals[force_key] = (
                        force_visual,
                        source_id,
                        interaction,
                        force_length,
                        original_color,
                        displayed_length,
                        individual_force,
                    )

            for source_id, force_visual in electron.individual_force_visuals.items():
                if source_id not in sources:
                    force_visual.visible = False

            if controls.electron_forces_visible and force_length > 0:
                selected_force_visuals[resultant_key] = (
                    electron.force_visual,
                    "resultierende Kraft",
                    "resultant",
                    force_length,
                    color.yellow,
                    displayed_length,
                    force,
                )

        for atom in atoms:
            force_direction = nuclear_resultants.get(
                f"{atom.atom_id}:nucleus",
                (0.0, 0.0, 0.0),
            )
            force_length = length(force_direction)
            arrow_length = min(0.85, force_length * 4.0)
            if force_length > 0:
                force_direction = scale(force_direction, arrow_length / force_length)
            atom.force_visual.pos = atom.visual.pos
            atom.force_visual.axis = vector(*force_direction)
            force_key = (
                f"{atom.atom_id}:nucleus",
                "resultant",
                "resultant",
            )
            atom.force_visual.color = (
                vector(0.56, 0.0, 1.0)
                if controls.selected_force_key == force_key and force_length > 0
                else color.yellow
            )
            if force_length > 0:
                selected_force_visuals[force_key] = (
                    atom.force_visual,
                    "resultierende Kraft",
                    "resultant",
                    force_length,
                    color.yellow,
                    arrow_length,
                    atom.applied_force,
                )

        if controls.selected_force_key not in selected_force_visuals:
            controls.selected_force_key = None

        for pair in lone_pairs:
            first, second = pair.electrons
            start = first.visual.pos
            end = second.visual.pos
            pair.visual.pos = start
            pair.visual.axis = end - start
