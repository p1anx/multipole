import numpy as np
import matplotlib.pyplot as plt

def rotate_vector(x, y, phi):
    # Rotate the vector (x, y) by phi
    x_rot = x * np.cos(phi) - y * np.sin(phi)
    y_rot = x * np.sin(phi) + y * np.cos(phi)
    return x_rot, y_rot

def convert_to_cartesian(r, theta):
    y = r * np.cos(theta)
    x = r * np.sin(theta)
    return x, y
def convert_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta
# convert cartesian coordinate x, y = ((2cos(theta), sin(theta))) to polar coordinate
def cartesian_to_polar_Vector(x, y, theta):
    length = len(x)
    cart2pol = np.zeros((length, 2))
    for i in range(length):
        vector = np.array([x[i], y[i]])
        cart2pol_vector = np.array([[np.sin(theta[i]), np.cos(theta[i])], [np.cos(theta[i]), -np.sin(theta[i])]])
        cart2pol[i] = np.dot(cart2pol_vector, vector)
    Br = cart2pol[:, 0]
    Btheta = cart2pol[:, 1]
    return Br, Btheta

def test_cartesian_to_polar_Vector():
    x, y = [[1, 0], [0,1]]
    theta = [0, np.pi]
    br, btheta = cartesian_to_polar_Vector(x, y, theta)
    print(br, btheta)

def test_cartesian_to_polar_Vector1():
    N = 100
    theta = np.linspace(np.pi/6, 2 * np.pi, N)
    r = np.linspace(0, 10, N)
    d = 0  # distance between poles
    theta, r = np.meshgrid(theta, r)
    x, y = convert_to_cartesian(r, theta)

    r1, theta1 = convert_to_polar(x-d, y)
    Bx1 = 3*np.cos(theta1)*np.sin(theta1)/r1**3   # Bx field
    By1 = (3*np.cos(theta1)**2-1)/r1**3   # By field
    bx1, by1 = cartesian_to_polar_Vector(Bx1, By1, theta1)
    plt.contourf(x, y, bx1)
    plt.show()

def polar():
    N = 10
    theta = np.linspace(np.pi/6, 2 * np.pi, N)
    r = np.linspace(0, 10, N)
    # r = 200  # radius
    d = 1  # distance between poles
    theta, r = np.meshgrid(theta, r)
    x, y = convert_to_cartesian(r, theta)

    r1, theta1 = convert_to_polar(x-d, y)
    # B1 = np.cos(theta1) / r1**3  # B field
    # B1 = 20*np.log10(abs(B1))  # Convert to dB
    # B = np.cos(theta) / r # B field

    Bx1 = 3*np.cos(theta1)*np.sin(theta1)/r1**3   # Bx field
    By1 = (3*np.cos(theta1)**2-1)/r1**3   # By field
    Bx1 = 20*np.log10(abs(Bx1))  # Convert to dB

    # print(Bx1)
    # B = np.cos(theta) / r # B field
    r2, theta2 = convert_to_polar(x+d, y)
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
    plt.axis('equal')
    plt.show()
if __name__ == "__main__":
    # polar()
    test_cartesian_to_polar_Vector1()
    # test_rotate_vector()
    # test_array_rotation()
