from math import cos, radians, sin, sqrt

Vec3 = tuple[float, float, float]


def add(first: Vec3, second: Vec3) -> Vec3:
    return tuple(first[index] + second[index] for index in range(3))


def subtract(first: Vec3, second: Vec3) -> Vec3:
    return tuple(first[index] - second[index] for index in range(3))


def scale(value: Vec3, scalar: float) -> Vec3:
    return tuple(component * scalar for component in value)


def dot(first: Vec3, second: Vec3) -> float:
    return sum(first[index] * second[index] for index in range(3))


def cross(first: Vec3, second: Vec3) -> Vec3:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def length(value: Vec3) -> float:
    return sqrt(dot(value, value))


def direction_towards(start: Vec3, end: Vec3, z_offset: float = 0.0) -> Vec3:
    dx, dy, dz = (end[index] - start[index] for index in range(3))
    dz += z_offset
    length = sqrt(dx * dx + dy * dy + dz * dz)
    return dx / length, dy / length, dz / length


def direction_at_angle(angle_degrees: float, z_offset: float) -> Vec3:
    angle = radians(angle_degrees)
    length = sqrt(1 + z_offset * z_offset)
    return cos(angle) / length, sin(angle) / length, z_offset / length


def remove_radial_velocity(velocity: Vec3, radial_direction: Vec3) -> Vec3:
    return subtract(velocity, scale(radial_direction, dot(velocity, radial_direction)))
