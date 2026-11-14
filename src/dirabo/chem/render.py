from collections import deque
from math import ceil

from vpython import color, curve, cylinder, rate, scene, sphere, vector

from dirabo.chem.keybd_handler import bind_keyboard_controls
from dirabo.chem.math import Vec3
from dirabo.chem.particles import (
    SIMULATION_DT,
    SIMULATION_STEPS_PER_FRAME,
    TRACE_MARKER_INTERVAL_SECONDS,
    TRACE_POINT_LIMIT,
    Molecule,
    MoleculeCluster,
    NuclearTrace,
    integrate_electrons,
    recenter_system,
)
from dirabo.chem.psystems import get_color


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
    show_molecules(MoleculeCluster(molecule))


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
        atom.visual = sphere(
            pos=vector(*atom.position),
            radius=atom.display_radius,
            color=vector(*get_color(atom.element)),
        )
        atom.normal_display_radius = atom.display_radius

    for bond in bonds:
        bond.visual = cylinder(
            pos=vector(0, 0, 0),
            axis=vector(0, 0, 0),
            radius=0.035,
            color=color.white,
        )

    for electron in electrons:
        electron.visual = sphere(
            pos=vector(*electron.position),
            radius=electron.display_radius,
            color=color.yellow,
        )

    for pair in lone_pairs:
        pair.visual = cylinder(
            pos=vector(0, 0, 0),
            axis=vector(0, 0, 0),
            radius=0.035,
            color=vector(0.82, 0.75, 0.62),
        )

    controls = bind_keyboard_controls(molecules, traces)
    trace_sample_index = 0
    while True:
        rate(60)
        if controls.paused:
            continue
        for _ in range(SIMULATION_STEPS_PER_FRAME):
            integrate_electrons(
                molecules,
                SIMULATION_DT,
                inter_molecular_bond_forces_enabled=controls.inter_molecular_bond_forces_enabled,
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

        for electron in electrons:
            electron.visual.pos = vector(*electron.position)

        for bond in bonds:
            first, second = bond.electrons
            start = first.visual.pos
            end = second.visual.pos
            bond.visual.pos = start
            bond.visual.axis = end - start

        for pair in lone_pairs:
            first, second = pair.electrons
            start = first.visual.pos
            end = second.visual.pos
            pair.visual.pos = start
            pair.visual.axis = end - start
