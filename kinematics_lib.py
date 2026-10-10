from dataclasses import dataclass

import sympy as sp

# set of dunctions defining elementar SO(3) transformations
def rotation_x(angle) -> sp.Matrix:
    return sp.Matrix([
        [1, 0, 0],
        [0, sp.cos(angle), -sp.sin(angle)],
        [0, sp.sin(angle), sp.cos(angle)]
    ])

def rotation_y(angle) -> sp.Matrix:
    return sp.Matrix([
        [sp.cos(angle), 0, sp.sin(angle)],
        [0, 1, 0],
        [-sp.sin(angle), 0, sp.cos(angle)]
    ])

def rotation_z(angle) -> sp.Matrix:
    return sp.Matrix([
        [sp.cos(angle), -sp.sin(angle), 0],
        [sp.sin(angle), sp.cos(angle), 0],
        [0, 0, 1]
    ])

def rotation_matrix(axis, angle) -> sp.Matrix:
    if axis == 'x':
        return rotation_x(angle)
    elif axis == 'y':
        return rotation_y(angle)
    elif axis == 'z':
        return rotation_z(angle)
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'.")

def rotation_matrix_to_se3(rotation_matrix) -> sp.Matrix:
    se3_matrix = sp.eye(4)
    se3_matrix[:3, :3] = rotation_matrix
    return se3_matrix

def rot(axis, angle) -> sp.Matrix:
    return rotation_matrix_to_se3(rotation_matrix(axis, angle))

def translation_matrix(translation_vector) -> sp.Matrix:
    se3_matrix = sp.eye(4)
    se3_matrix[:3, 3] = translation_vector
    return se3_matrix

trans = translation_matrix

@dataclass
class DenavitHartenbergFrameTransform:
    theta: sp.Symbol
    d: sp.Symbol
    l: sp.Symbol
    alpha: sp.Symbol

    # function to compute the transformation matrix using Denavit-Hartenberg parameters
    def transformation_matrix(self) -> sp.Matrix:
        return rot('z', self.theta) * trans(sp.Matrix([self.l, 0, self.d])) * rot('x', self.alpha)

    def latexify(self) -> str:
        return f"\\theta: {sp.latex(self.theta)}, d: {sp.latex(self.d)}, l: {sp.latex(self.l)}, \\alpha: {sp.latex(self.alpha)}"


def kinematics(transforms: list[DenavitHartenbergFrameTransform]) -> sp.Matrix:
    result = sp.eye(4)
    for transform in transforms:
        result *= transform.transformation_matrix()
    return result

@dataclass
class DenavitHartenbergWalkthroughStep:
    current_transform: DenavitHartenbergFrameTransform
    current_matrix: sp.Matrix
    total_matrix: sp.Matrix

def kinematics_walkthrough(transforms: list[DenavitHartenbergFrameTransform]) -> list[DenavitHartenbergWalkthroughStep]:
    result = sp.eye(4)
    walkthrough_steps = []
    for transform in transforms:
        current_matrix = transform.transformation_matrix()
        result *= current_matrix
        walkthrough_steps.append(DenavitHartenbergWalkthroughStep(
            current_transform=transform,
            current_matrix=current_matrix,
            total_matrix=result
        ))
    return walkthrough_steps

def kinematics_walkthrough_to_latex(walkthrough_steps: list[DenavitHartenbergWalkthroughStep]) -> str:
    latex_str = ""
    for step in walkthrough_steps:
        latex_str += f"Current Transform: {step.current_transform.latexify()}\n"
        latex_str += f"Current Matrix:\n{sp.latex(step.current_matrix)}\n"
        latex_str += f"Total Matrix:\n{sp.latex(step.total_matrix)}\n\n"
    return latex_str
