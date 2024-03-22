import numpy as np
import matplotlib.pyplot as plt


def rotate_vector(x, y, phi):
    # Rotate the vector (x, y) by phi
    x_rot = x * np.cos(phi) + y * np.sin(phi)
    y_rot = -x * np.sin(phi) + y * np.cos(phi)
    return x_rot, y_rot


def convert_to_cartesian(r, theta):
    y = r * np.cos(theta)
    x = r * np.sin(theta)
    return x, y


def convert_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(x, y)
    return r, theta


# convert cartesian coordinate x, y = ((2cos(theta), sin(theta))) to polar coordinate
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


def cartesian_to_polar_Vector1(x, y, theta):
    vector = np.array([x, y])
    cart2pol_vector = np.array(
        [[np.sin(theta), np.cos(theta)], [np.cos(theta), -np.sin(theta)]]
    )
    print(cart2pol_vector)
    cart2pol = np.dot(cart2pol_vector, vector.T)
    return cart2pol


def rotate_magnetic_dipole(r, theta, d, phi):
    x, y = convert_to_cartesian(r, theta)
    x1, y1 = rotate_vector(d, 0, phi)
    print(x1, y1)
    r1, theta1 = convert_to_polar(x - x1, y - y1)
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

    Bx2, By2 = rotate_vector(Bx1, By1, phi)
    Br, Btheta = cartesian_to_polar_Vector(Bx2, By2, theta)
    plt.plot(theta, Br)
    plt.show()


def rotate_magnetic_dipole1(r1, theta1):
    x, y = convert_to_cartesian(r1, theta1)
    # r1, theta1 = convert_to_polar(x, y)
    Bx11 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By11 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field
    # print('bx11:', Bx11)
    # print('by11:', By11)
    degree = np.pi / 5
    # degree = degree * 20
    Bx1, By1 = rotate_vector(Bx11, By11, degree)
    # print('bx1:', Bx1)
    # print('by1:', By1)
    Bx2, By2 = cartesian_to_polar_Vector(Bx1, By1, theta1 + degree)
    # print('bx2:', Bx2)
    # print('by2:', By2)
    # Bx2 = Bx1
    # By2 = By1

    ratio = (max(Bx2) - min(Bx2)) / (max(By2) - min(By2))
    maxBx2_value = max(Bx2)
    maxby2_value = max(By2)
    maxBx2_index = np.argmin(Bx2)
    maxBy2_index = np.argmin(By2)

    theta1 = np.degrees(theta1 + degree)
    print("degree:", np.degrees(degree))
    print("delta:", theta1[maxBx2_index] - theta1[maxBy2_index])
    plt.title("Bx/By = %f" % ratio)
    plt.plot(theta1, Bx2, label="Br")
    plt.plot(theta1, By2, label="Btheta")
    plt.legend()
    plt.show()


def init_multipole(r, theta, d, phi):
    x, y = convert_to_cartesian(r, theta)
    x1, y1 = rotate_vector(d, 0, phi)
    r1, theta1 = convert_to_polar(x - x1, y - y1)
    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field
    Bx2, By2 = rotate_vector(Bx1, By1, phi)
    Br, Btheta = cartesian_to_polar_Vector(Bx2, By2, theta1 + phi)
    plt.plot(theta1, Br, label="Br")
    plt.plot(theta1, Btheta, label="Btheta")
    plt.legend()
    plt.show()

    return Br, Btheta


def polar():
    N = 10
    theta = np.linspace(np.pi / 6, 2 * np.pi, N)
    r = np.linspace(0, 10, N)
    # r = 200  # radius
    d = 1  # distance between poles
    theta, r = np.meshgrid(theta, r)
    x, y = convert_to_cartesian(r, theta)

    r1, theta1 = convert_to_polar(x - d, y)
    # B1 = np.cos(theta1) / r1**3  # B field
    # B1 = 20*np.log10(abs(B1))  # Convert to dB
    # B = np.cos(theta) / r # B field

    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field
    Bx1 = 20 * np.log10(abs(Bx1))  # Convert to dB

    # print(Bx1)
    # B = np.cos(theta) / r # B field
    r2, theta2 = convert_to_polar(x + d, y)
    B2 = np.cos(theta2) / r2**3  # B field
    # B2 = 20*np.log10(abs(B2))  # Convert to dB
    # B = B1 - B2
    # B = 20*np.log10(abs(B))  # Convert to dB

    # plt.contourf(x, y, B1)
    plt.contourf(x, y, Bx1)
    # plt.contourf(x, y, B1, cmap='coolwarm')  # Use the 'coolwarm' colormap for better color differentiation

    # y = r * np.cos(theta)
    # x = r * np.sin(theta)
    # plt.plot(x, y, 'r-')
    plt.axis("equal")
    plt.show()


def test():
    theta = np.linspace(0, np.pi * 2, 5)
    r = 20  # radius
    x, y = convert_to_cartesian(r, theta)
    r1, theta1 = convert_to_polar(x, y)
    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    theta2 = np.degrees(theta1)

    print(theta2)
    a = np.linspace(0, 2, 5)
    print(a)
    plt.plot(theta2, a)
    plt.show()


if __name__ == "__main__":
    # polar()
    # test_cartesian_to_polar_Vector1()
    # test_rotate_vector()
    # test_array_rotation()
    test()
