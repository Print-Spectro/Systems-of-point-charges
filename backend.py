from numpy import *
from itertools import combinations
e0 = 8.854187817e-12
e = 1.6021765653e-19
r_x = lambda th: array(((1,0,0),(0,cos(th),-sin(th)),(0,sin(th),cos(th)))) #rotation matrices
r_y = lambda ph: array(((cos(ph),0, sin(ph)),(0,1,0),(-sin(ph),0,cos(ph))))
r_z = lambda az: array(((cos(az),-sin(az), 0),(sin(az),cos(az),0),(0,0,1)))


class charge:
    """
    Stores the position and charge. 
    """
    def  __init__(self, position, charge, **kwargs):
        if len(position) == 3:
           self.position = array(position)
        else:
            print("position is a 3 dimensional list")
        self.charge = charge
        self.colour = kwargs.get('colour', "b")

class charge_system:
    
    def __init__(self, charges):
        if isinstance([], type(charges)):    
            self.charges = charges
        else:
            print("charge_system object takes a list of charge objects")
        self.positions = [i.position for i in self.charges]
        self.position_factor = 1e-10
        self.charge_factor = 1.602176565e-19
    def potential(self, location, **kwargs):
        """
        
        Parameters
        ----------
        location : array or list
            3d coordinates at the position in space that you want to know the 
            potential of the system of charges
        location_factor : float
            the units for your unputted position coordinates. By defaut it's in
            angstroms

        Returns
        -------
        v : float
            the total potential at the inputted point in space

        """
        location_factor = kwargs.get('location_factor', 1e-10)
        location = array(location)
        v = 0
        for i in self.charges:
            if linalg.norm(location-i.position) == 0:
                print("Potential error; intersection with point charge")
                break
            else:
                v += i.charge / linalg.norm(location-i.position) 
        return v * self.charge_factor/(4*pi*e0*self.position_factor)
            
        
        ...
    def electric_field(self, location, **kwargs):
        """
        
        Parameters
        ----------
        location : TYPE 3 element array or list
            takes a 3d position coordinate 
        Location_factor : TYPE float
            the units of the inputted position coordinates in metres or whatever
            you want really

        Returns
        -------
        array
            unit vector for the net electric field
        float
            magnitude of the net electric field

        """
        location_factor = kwargs.get('location_factor', 1e-10)
        location = array(location)
        net = array([0.0,0.0,0.0])
        for i in self.charges:
            s = (location-i.position)
            radius = linalg.norm(s)
            if radius == 0:
                print("Location intersects with point charge; divide by zero")
                break
            net += i.charge*s/(radius**3)
        if linalg.norm(net) == 0.0:
           return (array([0,0,0]), 0)
        return net/linalg.norm(net), linalg.norm(net)*self.charge_factor/(4*pi*e0*self.position_factor**2)
        
        
        
    def potential_energy(self):
        """
        Returns
        -------
        float
            The potential electric energy of the system relative to moving each
            charge in from infinity.

        """
        combinator = lambda f, data: sum(list(map(f, combinations(data, 2)))) #this finds every pair of charges with no repeats, applys Coulomb's law and sums
        f = lambda a: a[0].charge*a[1].charge/(linalg.norm(a[0].position-a[1].position))
        #above is the function being applied to pairs of charges: The energy of every pair is calculated then summed in the combinator
        return combinator(f, self.charges)*e**2/(4*pi*e0*self.position_factor)
    
    def water(self, **kwargs):
        """
        Parameters
        ----------
        structure : TYPE list of charge objects
            This is where you define the structure of your molecule
        
        Optional parameters
        -------------------
        centre: TYPE 3 element array or list
            This translates the rotation centre in cartesian space; the structure
            is rotated about (0,0,0), then translated by
            the given value of the centre variable, which is by default 
            (0,0,0).
        radius : TYPE float
            This gives the distance of the inputted geometry from the centre
            of rotation. By default 0.
        x_angle, y_angle, z_angle : Type float 
            angle of rotation about each axis in radians. All 0 by default
        Returns
        -------
        Water at the position given by the parameters you input. By default,
        the water will lie with the oxygen at (0,0,0) with the hydrogens
        bisected by the positive x axis and lying in the x y plane.

        """
        th, ph, az = kwargs.get('x_angle', 0), kwargs.get('y_angle', 0), kwargs.get('z_angle', 0)
        centre, radius = kwargs.get('centre', [0,0,0]), kwargs.get('radius', 0)
        ox, h1, h2 = array([radius,0,0]), array([radius + 0.8565, 0.7527, 0]), array([radius + 0.8565, -0.7527, 0])
        for i in [r_x(th), r_y(ph), r_z(az)]:
            ox, h1, h2 = i.dot(ox), i.dot(h1), i.dot(h2)
           
        self.charges += [charge(ox + array(centre), -0.80628, colour = "r"), charge(h1 + array(centre), 0.40314, colour = "b"), charge(h2 + array(centre), 0.40314, colour = "b")]
    
    def molecule(self, structure, **kwargs):
        """
        Parameters
        ----------
        structure : TYPE list of charge objects
            This is where you define the structure of your molecule
        
        Optional parameters
        -------------------
        centre: TYPE 3 element array or list
            This translates the rotation centre in cartesian space; the structure
            is rotated about (0,0,0), then translated by
            the given value of the centre variable, which is by default 
            (0,0,0).
        radius : TYPE float
            This gives the distance of the inputted geometry from the centre
            of rotation. By default 0.
        x_angle, y_angle, z_angle : Type float 
            angle of rotation about each axis in radians. All 0 by default
        Returns
        -------
        The given structure rotated with the given radius about the centre will
        be appended to the charges variable so that calculations can be carried 
        out

        """
        th, ph, az = kwargs.get('x_angle', 0), kwargs.get('y_angle', 0), kwargs.get('z_angle', 0)
        centre, radius = kwargs.get('centre', [0,0,0]), kwargs.get('radius', 0)
        for i in structure:
            i.position[0] = i.position[0] + radius
            for j in [r_x(th), r_y(ph), r_z(az)]:
                i.position = j.dot(i.position)
        self.charges += structure
    def plot(self):
        """
        

        Returns
        -------
        Tuple to be plotted by matplotlib

        """
        return [i.position[0] for i in self.charges], [i.position[1] for i in self.charges], [i.position[2] for i in self.charges]
# class water:
#     def __init__(self, poistion, angles):
#         self.position = position
#         self.angles = angles
#         self.charges = [charge(self.position, )]

#box = charge_system([charge([2,0,0], 1), charge([-2,0,0], 1)])
# print(box.positions)

#water_ = [charge([0,0,0], -0.80628), charge([0.8565, 0.7527, 0], 0.40314), charge([0.8565, -0.7527, 0], 0.40314)]
#box.molecule(water_, radius = 1, z_angle = pi/2)
#box.water([0,0,0], 1, x_angle = pi/2, y_angle = 0)
#[print(i.position) for i in box.charges]

#print(box.electric_field([0,0,0]))  
#print("potential",box.potential([0,0,0], location_factor = 1e-10))
#print("energy",box.potential_energy())
# #If you don't pass `out` the indices where (b == 0) will be uninitialized!
# c = divide(a, b, out=zeros_like(a), where=b!=0)
# print(c)
# th = pi/2
# r_x = array(((1,0,0),(0,cos(th),-sin(th)),(0,sin(th),cos(th))))
# print(r_x.dot(array([0.8565, 0.7527, 0])))
