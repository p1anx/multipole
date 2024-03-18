
import numpy as np
import matplotlib.pyplot as plt

def rotate_vector(x, y, phi):
    # Rotate the vector (x, y) by phi
    x_rot = x * np.cos(phi) - y * np.sin(phi)
    y_rot = x * np.sin(phi) + y * np.cos(phi)
    return x_rot, y_rot
def test_array_rotation():
    x = np.array([1, 0])
    y = np.array([0, 1])
    phi = np.pi/2
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(x_rot, y_rot)  # [6.12323399e-17 1.00000000e+00] [ 1.00000000e+00 -1.22464680e-16]
    x = np.array([0, 1])
    y = np.array([1, 0])
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(x_rot, y_rot)  # [-1.  1.] [1.22464680e-16 6.12323399e-17]
def test_rotate_vector():
    x = 1
    y = 0
    phi = np.pi/2
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(x_rot, y_rot)  # 6.123233995736766e-17 1.0
    x = 0
    y = 1
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(x_rot, y_rot)  # -1.0 1.2246467991473532e-16
    if __name__ == "__main__":
        test_array_rotation()
