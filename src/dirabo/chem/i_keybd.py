"""Keyboard controls for the VPython chemistry visualization."""

from typing import Any

from vpython import scene, vector

from dirabo.chem.particles import (
    Atom,
    ELECTRON_DISPLAY_RADIUS,
    Electron,
    MoleculeCluster,
    NuclearTrace,
)

_SELECTED_PARTICLE_COLOR = vector(0.56, 0.0, 1.0)


class KeyboardControls:
    def __init__(
        self,
        molecules: MoleculeCluster,
        electrons_by_id: dict[str, Electron],
        traces: dict[str, NuclearTrace],
        nuclear_force_visuals: dict[tuple[str, ...], tuple[Any, Any]],
        selected_force_visuals: dict[
            tuple[str, ...],
            tuple[Any, str, str, float, Any, float, tuple[float, float, float]],
        ],
    ) -> None:
        self.molecules = molecules
        self.electrons_by_id = electrons_by_id
        self.traces = traces
        self.nuclear_force_visuals = nuclear_force_visuals
        self.selected_force_visuals = selected_force_visuals
        self.selectable_particles: tuple[tuple[str, Atom | Electron], ...] = tuple(
            (f"{atom.atom_id} (Kern)", atom) for atom in molecules.atoms
        ) + tuple((electron.electron_id, electron) for electron in electrons_by_id.values())
        self.selected_particle_index = -1
        self.selected_particle_original_color: Any = None
        self.selected_force_key: tuple[str, ...] | None = None
        self.electron_forces_visible = True
        self.nuclear_force_traces_visible = True
        self.small_nuclei = False
        self.traces_visible = True
        self.paused = False
        self.inter_molecular_bond_forces_enabled = True

    @property
    def selected_particle_id(self) -> str | None:
        if self.selected_particle_index < 0:
            return None
        _, particle = self.selectable_particles[self.selected_particle_index]
        if isinstance(particle, Atom):
            return f"{particle.atom_id}:nucleus"
        return particle.electron_id

    def handle_keydown(self, event: Any) -> None:
        key = str(event.key).lower()
        if key == "p":
            if not self.selectable_particles:
                print("Keine Partikel zur Auswahl")
                return

            self._restore_selected_force()
            if self.selected_particle_index >= 0:
                _, previous_particle = self.selectable_particles[self.selected_particle_index]
                previous_particle.visual.color = self.selected_particle_original_color

            self.selected_particle_index = (self.selected_particle_index + 1) % len(self.selectable_particles)
            name, particle = self.selectable_particles[self.selected_particle_index]
            current_color = particle.visual.color
            self.selected_particle_original_color = vector(
                current_color.x,
                current_color.y,
                current_color.z,
            )
            particle.visual.color = _SELECTED_PARTICLE_COLOR
            x, y, z = particle.position
            print(f"Partikel {name}: Position=({x:.3f}, {y:.3f}, {z:.3f})")
        elif key == "k":
            if self.selected_particle_index < 0:
                print("Zuerst ein Partikel mit 'p' auswählen")
                return
            selected_particle_id = self.selected_particle_id
            force_keys = [
                force_key for force_key in self.selected_force_visuals if force_key[0] == selected_particle_id
            ]
            if not force_keys:
                print("Für das ausgewählte Partikel sind keine Kraftpfeile visualisiert")
                return

            previous_force_key = self.selected_force_key
            self._restore_selected_force()
            if previous_force_key in force_keys:
                next_index = (force_keys.index(previous_force_key) + 1) % len(force_keys)
            else:
                next_index = 0
            self.selected_force_key = force_keys[next_index]
            (
                force_visual,
                source_id,
                interaction,
                magnitude,
                original_color,
                displayed_length,
                force_components,
            ) = self.selected_force_visuals[self.selected_force_key]
            force_visual.color = _SELECTED_PARTICLE_COLOR
            component_text = ", ".join(f"{component:.5g}" for component in force_components)
            if interaction == "resultant":
                print(
                    f"Resultierende Kraft an {selected_particle_id}: "
                    f"F=({component_text}), |F|={magnitude:.6g}, "
                    f"Pfeillänge={displayed_length:.3f}"
                )
            else:
                print(
                    f"Teil-Kraft an {selected_particle_id} von {source_id} "
                    f"({interaction}): F=({component_text}), |F|={magnitude:.6g}, "
                    f"Pfeillänge={displayed_length:.3f}"
                )
        elif key == "e":
            self.electron_forces_visible = not self.electron_forces_visible
            for electron in self.electrons_by_id.values():
                electron.force_visual.visible = self.electron_forces_visible
                for force_visual in electron.individual_force_visuals.values():
                    force_visual.visible = self.electron_forces_visible
            state = "an" if self.electron_forces_visible else "aus"
            print(f"Elektronenkräfte: {state}")
        elif key == "n":
            self.small_nuclei = not self.small_nuclei
            for atom in self.molecules.atoms:
                if atom.visual is not None:
                    atom.visual.radius = ELECTRON_DISPLAY_RADIUS if self.small_nuclei else atom.normal_display_radius
            state = "klein" if self.small_nuclei else "normal"
            print(f"Atomkerne: {state}")
        elif key == "t":
            self.traces_visible = not self.traces_visible
            for trace in self.traces.values():
                trace.visual.visible = self.traces_visible
                for marker_pair in trace.marker_visuals:
                    for marker in marker_pair:
                        marker.visible = self.traces_visible
            state = "an" if self.traces_visible else "aus"
            print(f"Atomkern-Traces: {state}")
        elif key == "a":
            self.nuclear_force_traces_visible = not self.nuclear_force_traces_visible
            for source_line, force_vector in self.nuclear_force_visuals.values():
                source_line.visible = self.nuclear_force_traces_visible
                force_vector.visible = self.nuclear_force_traces_visible
            state = "an" if self.nuclear_force_traces_visible else "aus"
            print(f"Nucleus-Kraft-Traces: {state}")
        elif key == "b":
            self.inter_molecular_bond_forces_enabled = not self.inter_molecular_bond_forces_enabled
            state = "an" if self.inter_molecular_bond_forces_enabled else "aus"
            print(f"Intermolekulare Kräfte gebundener Elektronen: {state}")
        elif key in (" ", "space"):
            self.paused = not self.paused
            state = "pausiert" if self.paused else "läuft weiter"
            print(f"Simulation {state}")

    def _restore_selected_force(self) -> None:
        if self.selected_force_key is None:
            return
        force_data = self.selected_force_visuals.get(self.selected_force_key)
        if force_data is not None:
            force_data[0].color = force_data[4]
        self.selected_force_key = None


def bind_keyboard_controls(
    molecules: MoleculeCluster,
    electrons_by_id: dict[str, Electron],
    traces: dict[str, NuclearTrace],
    nuclear_force_visuals: dict[tuple[str, ...], tuple[Any, Any]],
    selected_force_visuals: dict[
        tuple[str, ...],
        tuple[Any, str, str, float, Any, float, tuple[float, float, float]],
    ],
) -> KeyboardControls:
    controls = KeyboardControls(
        molecules,
        electrons_by_id,
        traces,
        nuclear_force_visuals,
        selected_force_visuals,
    )
    scene.bind("keydown", controls.handle_keydown)
    print(
        "Tastatur: p=Partikel auswählen, k=Kraftbeiträge durchlaufen, "
        "e=Elektronenkräfte, n=Atomgröße, t=Atomkern-Traces, "
        "a=Nucleus-Kraft-Traces, b=gebundene Elektron-Nucleus-Kräfte zwischen Molekülen, "
        "SPACE=Pause/Fortsetzen"
    )
    return controls
