"""Keyboard controls for the VPython chemistry visualization."""

from typing import Any

from vpython import scene

from dirabo.chem.particles import (
    ELECTRON_DISPLAY_RADIUS,
    MoleculeCluster,
    NuclearTrace,
)


class KeyboardControls:
    def __init__(
        self,
        molecules: MoleculeCluster,
        traces: dict[str, NuclearTrace],
    ) -> None:
        self.molecules = molecules
        self.traces = traces
        self.small_nuclei = False
        self.traces_visible = True
        self.paused = False
        self.inter_molecular_bond_forces_enabled = True

    def handle_keydown(self, event: Any) -> None:
        key = str(event.key).lower()
        if key == "n":
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
        elif key == "b":
            self.inter_molecular_bond_forces_enabled = not self.inter_molecular_bond_forces_enabled
            state = "an" if self.inter_molecular_bond_forces_enabled else "aus"
            print(f"Intermolekulare Kräfte gebundener Elektronen: {state}")
        elif key in (" ", "space"):
            self.paused = not self.paused
            state = "pausiert" if self.paused else "läuft weiter"
            print(f"Simulation {state}")


def bind_keyboard_controls(
    molecules: MoleculeCluster,
    traces: dict[str, NuclearTrace],
) -> KeyboardControls:
    controls = KeyboardControls(molecules, traces)
    scene.bind("keydown", controls.handle_keydown)
    print(
        "Tastatur: n=Atomgröße, t=Atomkern-Traces, "
        "b=gebundene Elektron-Nucleus-Kräfte zwischen Molekülen, SPACE=Pause/Fortsetzen"
    )
    return controls
