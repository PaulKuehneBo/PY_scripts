import json
import numpy as np
from collections import defaultdict, deque
from dataclasses import dataclass

# ============================================================
# GROUP LIBRARY
# ============================================================

GROUP_LIBRARY = {
    "methyl": {
        "atoms": {"C": "C", "H1": "H", "H2": "H", "H3": "H"},
        "ports": {"carbon": "C"},
        "internal_bonds": [("C", "H1"), ("C", "H2"), ("C", "H3")],
    },
    "oxygen": {"atoms": {"O": "O"}, "ports": {"left": "O", "right": "O"}, "internal_bonds": []},
}


def update_group_library(new_groups):
    GROUP_LIBRARY.update(new_groups)


# ============================================================
# AtomPlacement
# ============================================================
@dataclass
class AtomPlacement:
    id: str
    element: str
    position: np.ndarray


# ============================================================
# CHEMICAL BOND LENGTHS
# ============================================================

BOND_LENGTH = {
    ("C", "H"): 1.09,
    ("H", "C"): 1.09,
    ("C", "O"): 1.43,
    ("O", "C"): 1.43,
    ("C", "C"): 1.54,
    ("O", "H"): 0.96,
    ("H", "O"): 0.96,
}


def _get_bond_length(a, b):

    return BOND_LENGTH.get((a, b), 1.5)


# ============================================================
# GEOMETRY TEMPLATES
# ============================================================

LINEAR = np.array(
    [
        [1, 0, 0],
        [-1, 0, 0],
    ],
    dtype=float,
)

TRIGONAL = np.array(
    [
        [1, 0, 0],
        [-0.5, 0.866, 0],
        [-0.5, -0.866, 0],
    ],
    dtype=float,
)

TETRA = np.array(
    [
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
    ],
    dtype=float,
)

for i in range(len(TETRA)):
    TETRA[i] /= np.linalg.norm(TETRA[i])

# ============================================================
# GEOMETRY SELECTION
# ============================================================


def _geometry(child_count):

    if child_count <= 1:
        return LINEAR[:1]

    if child_count == 2:
        return LINEAR

    if child_count == 3:
        return TRIGONAL

    return TETRA[:child_count]


# ============================================================
# VECTOR HELPERS
# ============================================================


def _normalize(vec: np.ndarray) -> np.ndarray:

    vec_length = np.linalg.norm(vec)

    if vec_length < 1e-8:
        return vec

    return vec / vec_length


def _perpendicular(vec: np.ndarray) -> np.ndarray:

    if abs(vec[0]) < 0.9:
        ref = np.array([1, 0, 0], dtype=float)
    else:
        ref = np.array([0, 1, 0], dtype=float)

    return _normalize(np.cross(vec, ref))


def _rotation_matrix_from_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:

    a = _normalize(a)
    b = _normalize(b)

    cross = np.cross(a, b)
    dot = np.dot(a, b)

    if dot > 0.99999:
        return np.eye(3)

    if dot < -0.99999:
        axis = _perpendicular(a)
        x, y, z = axis
        return np.array(
            [
                [2 * x * x - 1, 2 * x * y, 2 * x * z],
                [2 * x * y, 2 * y * y - 1, 2 * y * z],
                [2 * x * z, 2 * y * z, 2 * z * z - 1],
            ]
        )

    s = np.linalg.norm(cross)

    vx = np.array(
        [
            [0, -cross[2], cross[1]],
            [cross[2], 0, -cross[0]],
            [-cross[1], cross[0], 0],
        ]
    )

    return np.eye(3) + vx + vx @ vx * ((1 - dot) / (s * s))


# ============================================================
# JSON -> ATOM GRAPH
# ============================================================


def _build_atom_graph(molecule):

    atom_types = {}

    bonds = []

    port_map = {}

    for component in molecule["components"]:

        component_id = component["id"]
        group_name = component["group"]

        template = GROUP_LIBRARY[group_name]

        # create global atoms for that component
        for template_types in template["atoms"]:
            atom_id = template_types["id"]
            element = template_types["element"]

            global_atom = f"{component_id}:{atom_id}"
            atom_types[global_atom] = element

        # create ports
        for port_name, atom_name in template["ports"].items():

            port_map[(component_id, port_name)] = f"{component_id}:{atom_name}"

        # build internal bonds
        bond_list = template.get("bonds", [])
        for bond in bond_list:
            from_atom = bond["from"]
            to_atom = bond["to"]
            bonds.append((f"{component_id}:{from_atom}", f"{component_id}:{to_atom}"))

    # create inter-component bonds
    for connection in molecule["connections"]:

        from_component = connection["from"]["component"]
        from_port = connection["from"]["port"]

        to_component = connection["to"]["component"]
        to_port = connection["to"]["port"]

        a = port_map[(from_component, from_port)]
        b = port_map[(to_component, to_port)]

        bonds.append((a, b))

    return atom_types, bonds


# ============================================================
# GRAPH
# ============================================================


def _build_graph(bonds):

    graph = defaultdict(list)

    for a, b in bonds:

        graph[a].append(b)
        graph[b].append(a)

    return dict(graph)


def _build_tree(graph, root):

    parent = {root: None}
    children = defaultdict(list)

    queue = deque([root])

    while queue:

        current = queue.popleft()

        for neighbour in graph[current]:

            if neighbour in parent:
                continue

            parent[neighbour] = current
            children[current].append(neighbour)

            queue.append(neighbour)

    return dict(children)


# ============================================================
# LAYOUT
# ============================================================


def _place_recursive(atom_id, atom_types, children, positions, direction):

    childs = children.get(atom_id, [])

    if not childs:
        return

    local_vectors = _geometry(len(childs))

    rotation = _rotation_matrix_from_vectors(np.array([0, 0, 1]), _normalize(direction))
    rotated = [rotation @ vec for vec in local_vectors]

    for child, vec in zip(childs, rotated):
        length = _get_bond_length(atom_types[atom_id], atom_types[child])
        positions[child] = positions[atom_id] + _normalize(vec) * length
        _place_recursive(child, atom_types, children, positions, vec)


def _generate_atom_positions(atom_types, bonds):

    graph = _build_graph(bonds)

    root = max(graph.keys(), key=lambda n: len(graph[n]))

    children = _build_tree(graph, root)

    positions = {root: np.array([0.0, 0.0, 0.0])}

    _place_recursive(root, atom_types, children, positions, np.array([0.0, 0.0, 1.0]))

    return positions


# ============================================================
# PUBLIC API
# ============================================================


def layout_molecule(molecule) -> list[AtomPlacement]:

    atom_types, bonds = _build_atom_graph(molecule)

    positions = _generate_atom_positions(atom_types, bonds)

    result = []
    for atom_id, position in positions.items():
        atom_placement = AtomPlacement(id=atom_id, element=atom_types[atom_id], position=position)
        result.append(atom_placement)

    return result


def load_layout(json_path: str, verbose: bool = False) -> dict[str, list[AtomPlacement]]:
    with open(json_path, "r", encoding="utf-8") as f:
        soup = json.load(f)

    update_group_library(soup["groups"])

    bom = {}
    molecules = soup["molecules"]
    for name, molecule in molecules.items():

        atoms: list[AtomPlacement] = layout_molecule(molecule)
        if verbose:
            print()
            print("MOLECULE:", name)

            for atom in atoms:
                print(atom.id, atom.element, np.round(atom.position, 3))

        bom[name] = atoms

    return bom


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    bom = load_layout("D:/SB/My Python Scripts/devl/molecular/molecule_simple_soup.json")

    print("\nDONE")
