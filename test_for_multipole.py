import multipole as mp
import numpy as np
import matplotlib.pyplot as plt


# def rotate_vector(x, y, phi):
#     # Rotate the vector (x, y) by phi
#     x_rot = x * np.cos(phi) - y * np.sin(phi)
#     y_rot = x * np.sin(phi) + y * np.cos(phi)
#     return x_rot, y_rot
def convert_to_cartesian(r, theta):
    y = r * np.cos(theta)
    x = r * np.sin(theta)
    return x, y


def convert_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta


# *********************************************************
def test_array_rotation():
    x = np.array([1, 0])
    y = np.array([0, 1])
    phi = np.pi / 2
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(
        x_rot, y_rot
    )  # [6.12323399e-17 1.00000000e+00] [ 1.00000000e+00 -1.22464680e-16]
    x = np.array([0, 1])
    y = np.array([1, 0])
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(x_rot, y_rot)  # [-1.  1.] [1.22464680e-16 6.12323399e-17]


def test_rotate_vector():
    x = 1
    y = 0
    phi = np.pi / 2
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(x_rot, y_rot)  # 6.123233995736766e-17 1.0
    x = 0
    y = 1
    x_rot, y_rot = rotate_vector(x, y, phi)
    print(x_rot, y_rot)  # -1.0 1.2246467991473532e-16


# *********************************************************
def cartesian_to_polar_Vector(x, y, theta):
    length = len(x)
    cart2pol = np.zeros((length, 2))
    for i in range(length):
        vector = np.array([x[i], y[i]])
        cart2pol_vector = np.array(
            [
                [np.sin(theta[i]), np.cos(theta[i])],
                [np.cos(theta[i]), -np.sin(theta[i])],
            ]
        )
        cart2pol[i] = np.dot(cart2pol_vector, vector)
    Br = cart2pol[:, 0]
    Btheta = cart2pol[:, 1]
    return Br, Btheta


def test_cartesian_to_polar_Vector():
    x, y = [[1, 0], [0, 1]]
    theta = [0, np.pi]
    br, btheta = cartesian_to_polar_Vector(x, y, theta)
    print(br, btheta)


def test_cartesian_to_polar_Vector1():
    N = 100
    theta = np.linspace(np.pi / 6, 2 * np.pi, N)
    r = np.linspace(0, 10, N)
    d = 0  # distance between poles
    theta, r = np.meshgrid(theta, r)
    x, y = convert_to_cartesian(r, theta)

    r1, theta1 = convert_to_polar(x - d, y)
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field
    bx1, by1 = cartesian_to_polar_Vector(Bx1, By1, theta1)
    plt.contourf(x, y, bx1)
    plt.show()


# ***********************************************************
def rotate_vector(x, y, phi):
    # Rotate the vector (x, y) by phi
    x_rot = x * np.cos(phi) + y * np.sin(phi)
    y_rot = -x * np.sin(phi) + y * np.cos(phi)
    return x_rot, y_rot


def rotate_magnetic_dipole(x, y, r, theta, phi):
    Bx1 = 3 * np.cos(theta) * np.sin(theta) / r**3  # Bx field
    By1 = (3 * np.cos(theta) ** 2 - 1) / r**3  # By field
    print("before:", Bx1, By1)
    Bx2, By2 = mp.rotate_vector(Bx1, By1, phi)
    print("after:", Bx2, By2)


def test_rotate_magnetic_dipole():
    theta = np.linspace(0, 2 * np.pi, 100)
    r = 20
    d = 0
    phi = 0
    mp.rotate_magnetic_dipole(r, theta, d, phi)


def test_rotate_magnetic_dipole1():
    theta = np.linspace(0, 2 * np.pi, 100)
    # print(np.degrees(theta))
    r = 1
    mp.rotate_magnetic_dipole1(r, theta)


def test_convert_vector():
    x, y = mp.rotate_vector(0, -1, np.pi / 2)
    print(x, y)
    x = 0
    y = -1
    a, b = mp.cartesian_to_polar_Vector1(x, y, np.pi / 2)
    # print(x, y)
    print(a, b)


def test_convert():
    theta = np.linspace(0, 2 * np.pi, 10)
    r = 20
    print("theta:", theta)
    x, y = mp.convert_to_cartesian(r, theta)
    print("x:", x)
    print("y:", y)
    r, theta = mp.convert_to_polar(x, y)
    print("r:", r)
    print("theta1:", theta)


def test_init_multipole():
    theta = np.linspace(0, 2 * np.pi, 100)
    r = 20
    d = 1
    phi3 = np.pi / 3
    Br3, Btheta3 = mp.init_multipole(r, theta, d, phi3)
    # plt.plot(theta, Br, label = 'Br')
    # plt.plot(theta, Btheta, label = 'Btheta')
    # plt.legend()
    # plt.show()


if __name__ == "__main__":
    # test_convert()
    # test_rotate_magnetic_dipole1()
    test_init_multipole()
