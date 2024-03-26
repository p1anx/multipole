import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import argrelextrema


path1 = '/home/clockx/Documents/code/multipole/figure/'

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

def test1():
    theta = np.linspace(0, np.pi * 2, 100)
    r = 20
    d = 5
    length = len(theta)
    x, y = convert_to_cartesian(r, theta)
    r1, theta1 = convert_to_polar(x-d, y)

    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

    Bx2, By2 = rotate_vector(Bx1, By1, np.pi/2)
    Bx = Bx1 + Bx2
    By = By1 + By2
    Br, Btheta = cartesian_to_polar_Vector(Bx, By, theta)
    # a = np.linspace(0, length-1, length)
    a = np.degrees(theta)
    plt.figure(1)
    plt.plot(a, Bx1, label = 'bx1')
    plt.plot(a, By1, label = 'by1')
    plt.legend()
    # plt.plot(a, np.degrees(theta1))
    # plt.axis("equal")
    plt.figure(2)
    plt.plot(a, Bx2, label = 'bx2')
    plt.plot(a, By2, label = 'by2')
    plt.legend()

    plt.figure(3)
    plt.plot(a, Bx, label = 'bx')
    plt.plot(a, By, label = 'by')
    plt.legend()
    plt.show()

def test2():
    theta = np.linspace(0, np.pi * 2, 100)
    r = 20
    d = 5
    length = len(theta)
    x, y = convert_to_cartesian(r, theta)
    r1, theta1 = convert_to_polar(x-d, y)

    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

    r2, theta2 = convert_to_polar(x+d, y)
    Bx2 = -3*np.cos(theta2)*np.sin(theta2)/r2**3
    By2 = -(3*np.cos(theta2)**2 - 1)/r2**3

    Bx = Bx1 + Bx2
    By = By1 + By2

    Br, Btheta = cartesian_to_polar_Vector(Bx, By, theta)
    a = np.degrees(theta)
    plt.figure(1)
    plt.plot(a, Bx1, label = 'Bx1')
    plt.plot(a, By1, label = 'By1')
    plt.legend()

    plt.figure(2)
    plt.plot(a, Bx2, label = 'Bx2')
    plt.plot(a, By2, label = 'By2')
    plt.legend()

    plt.figure(3)
    plt.plot(a, Bx, label = 'Bx')
    plt.plot(a, By, label = 'By')
    plt.legend()

    plt.figure(4)
    plt.plot(a, Br, label = 'Br')
    plt.plot(a, Btheta, label = 'Btheta')
    plt.legend()

    plt.show()
    
def test3():
    N = 100
    theta = np.linspace(0, np.pi * 2*(1-1/N), N)

    r = 20
    d = 5
    rot_angle = np.pi
    length = len(theta)
    x, y = convert_to_cartesian(r, theta)
    r1, theta1 = convert_to_polar(x-d, y)

    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

    Bx2, By2 = rotate_vector(Bx1, By1,  rot_angle)
    mid_index = round(length *(rot_angle/(2*np.pi)))+1
    
    # aa = Bx2[mid_index:] + Bx2[:mid_index]
    Bx2 = np.concatenate((Bx2[mid_index:], Bx2[:mid_index]))
    By2 = np.concatenate((By2[mid_index:], By2[:mid_index]))

    Bx = Bx1 + Bx2
    By = By1 + By2
    Br, Btheta = cartesian_to_polar_Vector(Bx, By, theta)
    # a = np.linspace(0, length-1, length)
    a = np.degrees(theta)
    plt.figure(1)
    plt.plot(a, Bx1, label = 'bx1')
    plt.plot(a, By1, label = 'by1')
    plt.legend()
    # plt.plot(a, np.degrees(theta1))
    # plt.axis("equal")
    plt.figure(2)
    plt.plot(a, Bx2, label = 'bx2')
    plt.plot(a, By2, label = 'by2')
    plt.legend()

    plt.figure(3)
    plt.plot(a, Bx, label = 'bx')
    plt.plot(a, By, label = 'by')
    plt.legend()

    plt.figure(4)
    plt.plot(a,Br, label = 'Br')
    plt.plot(a,Btheta, label = 'Btheta')
    plt.legend()

    path1 = '/home/clockx/Documents/code/multipole/figure/'
    plt.savefig(path1 +'test3.png')
    plt.show()

def test4():
# path1 = '/home/clockx/Documents/code/multipole/figure/'
    N = 100
    theta = np.linspace(0, np.pi * 2*(1-1/N), N)

    r = 20
    d = 5
    rot_angle = np.pi/2
    length = len(theta)
    x, y = convert_to_cartesian(r, theta)
    r1, theta1 = convert_to_polar(x-d, y)

    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field


    #rotate 90 degree
    Bx2, By2 = rotate_vector(Bx1, By1,  rot_angle)
    mid_index = length - round(length *(rot_angle/(2*np.pi)))
    
    # aa = Bx2[mid_index:] + Bx2[:mid_index]
    Bx2 = np.concatenate((Bx2[mid_index:], Bx2[:mid_index]))
    By2 = np.concatenate((By2[mid_index:], By2[:mid_index]))
    

    
    #rotatte 180 degree
    Bx3, By3 = rotate_vector(Bx1, By1, 2*rot_angle)
    mid_index = length - round(length *(2*rot_angle/(2*np.pi)))
    Bx3 = np.concatenate((Bx3[mid_index:], Bx3[:mid_index]))
    By3 = np.concatenate((By3[mid_index:], By3[:mid_index]))

    #rotatte 270 degree
    Bx4, By4 = rotate_vector(Bx1, By1, 3*rot_angle)
    mid_index = length - round(length *(3*rot_angle/(2*np.pi)))
    Bx4 = np.concatenate((Bx4[mid_index:], Bx4[:mid_index]))
    By4 = np.concatenate((By4[mid_index:], By4[:mid_index]))


    # Bx = Bx2 + Bx4
    # By = By2 + By4
    Bx = Bx1 + Bx2 + Bx3 + Bx4
    By = By1 + By2 + By3 + By4
    # Bx1, By1 = cartesian_to_polar_Vector(Bx1, By1, theta)
    # Bx2, By2 = cartesian_to_polar_Vector(Bx2, By2, theta)
    # Bx3, By3 = cartesian_to_polar_Vector(Bx3, By3, theta)
    # Bx4, By4 = cartesian_to_polar_Vector(Bx4, By4, theta)
    Br, Btheta = cartesian_to_polar_Vector(Bx, By, theta)
    # a = np.linspace(0, length-1, length)
    a = np.degrees(theta)
    plt.figure(1)
    plt.plot(a, Bx1, label = 'Bx1')
    plt.plot(a, By1, label = 'By1')
    plt.legend()
    plt.savefig(path1 + 'B1.png')

    # plt.plot(a, np.degrees(theta1))
    # plt.axis("equal")
    plt.figure(2)
    plt.plot(a, Bx2, label = 'Bx2')
    plt.plot(a, By2, label = 'By2')
    plt.legend()
    plt.savefig(path1 + 'B2.png')

    plt.figure(5)
    plt.plot(a, Bx3, label = 'Bx3')
    plt.plot(a, By3, label = 'By3')
    plt.legend()
    plt.savefig(path1 + 'B3.png')

    plt.figure(6)
    plt.plot(a, Bx4, label = 'Bx4')
    plt.plot(a, By4, label = 'By4')
    plt.legend()
    plt.savefig(path1 + 'B4.png')

    plt.figure(3)
    plt.plot(a, Bx, label = 'Bx')
    plt.plot(a, By, label = 'By')
    plt.legend()
    plt.savefig(path1 + 'Bx.png')

    plt.figure(4)
    plt.plot(a,Br, label = 'Br')
    plt.plot(a,Btheta, label = 'Btheta')
    plt.legend()
    plt.savefig(path1 + 'Br_Btheta.png')

    # plt.show()

def test5():
    N = 100
    theta = np.linspace(0, np.pi * 2*(1-1/N), N)

    r = 20
    d = 5
    Nm = 2 # number of magnet
    rot_angle = 2*np.pi/Nm
    length = len(theta)
    x, y = convert_to_cartesian(r, theta)

    for i in range(Nm):
        x1, y1 = rotate_vector(d, 0, i*rot_angle)
        r1, theta1 = convert_to_polar(x-x1, y-y1)

        for i in range(len(theta1)):
            if theta1[i] < 0:
                theta1[i] = theta1[i] + np.pi * 2
    
        Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
        By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

        

        Bx2 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
        By2 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

    path1 = '/home/clockx/Documents/code/multipole/figure/'
    plt.savefig(path1 + 'test3.png')
    plt.show()

def test6():
    theta = np.linspace(0, np.pi * 2, 100)
    r = 20
    d = 5
    length = len(theta)
    x, y = convert_to_cartesian(r, theta)
    r1, theta1 = convert_to_polar(x-d, y)

    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    
    Bx1 = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By1 = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

    #rotate 180
    r2, theta2 = convert_to_polar(x+d, y)
    Bx2 = -3*np.cos(theta2)*np.sin(theta2)/r2**3
    By2 = -(3*np.cos(theta2)**2 - 1)/r2**3

    #rotate 90
    r3, theta3 = convert_to_polar(x, y+d)
    By3 = -3*np.cos(theta3)*np.sin(theta3)/r3**3
    Bx3 = (3*np.cos(theta3)**2 - 1)/r3**3

    #rotate 270
    r4, theta4 = convert_to_polar(x, y-d)
    By4 = 3*np.cos(theta4)*np.sin(theta4)/r4**3
    Bx4 = -(3*np.cos(theta4)**2 - 1)/r4**3

    Bx = Bx1 + Bx2 + Bx3 + Bx4
    By = By1 + By2 + By3 + By4

    Br, Btheta = cartesian_to_polar_Vector(Bx, By, theta)
    a = np.degrees(theta)
    plt.figure(1)
    plt.plot(a, Bx1, label = 'Bx1')
    plt.plot(a, By1, label = 'By1')
    plt.legend()
    plt.savefig(path1 + 'B1.png')

    
    plt.figure(5)
    plt.plot(a, Bx3, label = 'Bx2')
    plt.plot(a, By3, label = 'By2')
    plt.legend()
    plt.savefig(path1 + 'B2.png')

    plt.figure(2)
    plt.plot(a, Bx2, label = 'Bx3')
    plt.plot(a, By2, label = 'By3')
    plt.legend()
    plt.savefig(path1 + 'B3.png')

    plt.figure(6)
    plt.plot(a, Bx4, label = 'Bx4')
    plt.plot(a, By4, label = 'By4')
    plt.legend()
    plt.savefig(path1 + 'B4.png')

    plt.figure(3)
    plt.plot(a, Bx, label = 'Bx')
    plt.plot(a, By, label = 'By')
    plt.legend()

    plt.figure(4)
    plt.plot(a, Br, label = 'Br')
    plt.plot(a, Btheta, label = 'Btheta')
    plt.legend()

    plt.show()
    
#Nm is number of magnet, d is distance from magnet to center point, 
#r is radius of magnetic probe
def rot_multipole(Nm, d, r, path): 
    if not os.path.exists(path):
        os.makedirs(path)

    fig_name = 'N=' + str(Nm) + ' d=' + str(d)+' r='+ str(r) + '.png'
    csv_name = 'N=' + str(Nm) + ' d=' + str(d)+' r='+ str(r) + '.csv'

# path1 = '/home/clockx/Documents/code/multipole/figure/'
    N = 100
    theta = np.linspace(0, np.pi * 2*(1-1/N), N)

    # r = 100
    rot_angle = np.pi*2/Nm
    length = len(theta)
    x, y = convert_to_cartesian(r, theta)
    r1, theta1 = convert_to_polar(x-d, y)

    for i in range(len(theta1)):
        if theta1[i] < 0:
            theta1[i] = theta1[i] + np.pi * 2
    
    Bx = 3 * np.cos(theta1) * np.sin(theta1) / r1**3  # Bx field
    By = (3 * np.cos(theta1) ** 2 - 1) / r1**3  # By field

    a = np.degrees(theta)
    sum_magnets = np.zeros((2, N))
    for i in range(Nm):
        Bx1, By1 = rotate_vector(Bx, By, i*rot_angle)
        mid_index = length - round(length *(i*rot_angle/(2*np.pi)))
        Bx1 = np.concatenate((Bx1[mid_index:], Bx1[:mid_index]))
        By1 = np.concatenate((By1[mid_index:], By1[:mid_index]))
        sum_magnets[0,:] = sum_magnets[0,:] + Bx1
        sum_magnets[1,:] = sum_magnets[1,:] + By1

        Bxx = sum_magnets[0, :]
        Byy = sum_magnets[1, :]

        Br, Btheta = cartesian_to_polar_Vector(Bxx, Byy, theta)
        Bdata = {'Br':Br, 'Btheta':Btheta}
        df = pd.DataFrame(Bdata)
        df.to_csv(path + csv_name, index = False)

        plt.figure(i)
        plt.plot(a,Br, label = 'Br')
        plt.plot(a,Btheta, label = 'Btheta')
        plt.xlabel('x/degree')
        plt.ylabel('y/B')
        plt.legend()
        # plt.savefig(path1 + 'Br_Btheta.png')

    # ratio = (max(Br) - min(Br)) / (max(Btheta) - min(Btheta))

    ratio = max(Br) / max(Btheta)
    # Br_max_index = argrelextrema(Br, np.greater)
    # Br_min_index = argrelextrema(Br, np.less)
    # Btheta_max_index = argrelextrema(Btheta, np.greater)
    # Btheta_min_index = argrelextrema(Btheta, np.less)
    #
    # print('brmaxindex:', Br_max_index)
    # br_max = np.mean(Br[Br_max_index])
    # print(br_max)
    # br_min = np.mean(Br[Br_min_index])
    # print(br_min)
    # btheta_max = np.mean(Btheta[Btheta_max_index])
    # btheta_min = np.mean(Btheta[Btheta_min_index])
    # ratio = (br_max - br_min) / (btheta_max - btheta_min)
    
    print('N='+str(Nm)+' ratio=', ratio)
    plt.savefig(path + fig_name)

    #rotate 90 degree
    # plt.show()

def experiment_multipole_N_influence():
    path = '/home/clockx/Documents/code/multipole/figure_r_less_d_ratio0.01/'
    N = 9
    ratio = 0.01
    d = 100
    r = d * ratio
    for i in range(1,N):
        rot_multipole(i, d, r, path)


    

def test():
    
    path = '/home/clockx/Documents/code/multipole/figure_ratio100_influence/'
    file_name = 'N=1 d=1 r=100.csv'
    data = pd.read_csv(path + file_name)
    Br = data['Br']
    Btheta = data['Btheta']
    plt.plot(range(len(Br)), Br)
    plt.plot(range(len(Br)), Btheta)
    plt.show()
    # print(data['Br'])

if __name__ == "__main__":
    # polar()
    # test_cartesian_to_polar_Vector1()
    # test_rotate_vector()
    # test_array_rotation()
    # test()
    experiment_multipole_N_influence()

