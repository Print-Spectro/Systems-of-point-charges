import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
e0 = 8.854187817e-12
e = 1.6021765653e-19
r_x = lambda th: np.array(((1,0,0),(0,np.cos(th),-np.sin(th)),(0,np.sin(th),np.cos(th)))) #x axis rotation matrix
r_y = lambda ph: np.array(((np.cos(ph),0, np.sin(ph)),(0,1,0),(-np.sin(ph),0,np.cos(ph))))
r_z = lambda az: np.array(((np.cos(az),-np.sin(az), 0),(np.sin(az),np.cos(az),0),(0,0,1)))

class charge:
    """
    Stores the position and charge. 
    """
    def  __init__(self, position, charge, **kwargs):
        if len(position) == 3:
           self.position = np.array(position)
        else:
            print("position is a 3 dimensional list")
        self.charge = charge
        self.colour = kwargs.get('colour', "black")
          
class charge_system:   
    def __init__(self, charges, **kwargs):
        """
        charge_system object stores charges and can do lots of things with them

        Parameters
        ----------
        charges : List of charge objects
            
        position_factor : float
            This gives the units of the coordinates, used for caluclating in si
            units
        charge_factor : float
            this is the units in C of the given charge, used for doing calculations
            in si units.
        Returns
        -------
        None.

        """
        if isinstance([], type(charges)):    
            self.charges = charges
        else:
            print("charge_system object takes a list of charge objects")
        #self.positions = [i.position for i in self.charges]
        self.position_factor = kwargs.get("position_factor", 1e-10)
        self.charge_factor = kwargs.get("charge_factor", 1.602176565e-19)
        self.bonds = []
    def potential(self, location, **kwargs):
        """
        
        Parameters
        ----------
        location : array or list
            3d coordinates at the position in space that you want to know the 
            potential of the system of charges

        Returns
        -------
        v : float
            the total potential at the inputted point in space

        """
        location = np.array(location)
        v = 0
        for i in self.charges:
            if np.linalg.norm(location-i.position) == 0:
                print("Potential error; intersection with point charge")
                break
            else:
                v += i.charge / np.linalg.norm(location-i.position) 
        return v * self.charge_factor/(4*np.pi*e0*self.position_factor)
            
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
        
        location = np.array(location)
        net = np.array([0.0,0.0,0.0])
        for i in self.charges:
            s = (location-i.position)
            radius = np.linalg.norm(s)
            if radius == 0:
                print("Location intersects with point charge; divide by zero")
                break
            net += i.charge*s/(radius**3)
        if np.linalg.norm(net) == 0.0:
           return (np.array([0,0,0]), 0)
        return net/np.linalg.norm(net), np.linalg.norm(net)*self.charge_factor/(4*np.pi*e0*self.position_factor**2)
            
    def potential_energy(self):
        """
        Returns
        -------
        float
            The potential electric energy of the system relative to moving each
            charge in from infinity.

        """
        combinator = lambda f, data: sum(list(map(f, combinations(data, 2)))) #this finds every pair of charges with no repeats, applys Coulomb's law and sums
        f = lambda a: a[0].charge*a[1].charge/(np.linalg.norm(a[0].position-a[1].position))
        #above is the function being applied to pairs of charges: The energy of every pair is calculated then summed in the combinator
        return combinator(f, self.charges)*e**2/(4*np.pi*e0*self.position_factor)
    
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
        centre, radius = np.array(kwargs.get('centre', [0,0,0])), kwargs.get('radius', 0)
        ox, h1, h2 = np.array([radius,0,0]), np.array([radius + 0.5865, 0.7572, 0]), np.array([radius + 0.5865, -0.7572, 0])

        for i in [r_x(th), r_y(ph), r_z(az)]:
            ox, h1, h2 = i.dot(ox), i.dot(h1), i.dot(h2)
           
        self.charges += [charge(ox  + centre, -0.658, colour = "r"), charge(h1 + centre , 0.329, colour = "b"), charge(h2  + centre, 0.40314, colour = "b")]
        self.bonds += [(ox + centre, h1 + centre), (ox + centre, h2 + centre)]
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
                i.position = j.dot(i.position) + np.array(centre)
        self.charges += structure

    def get_bonds(self):
        """
        
        Returns
        -------
        list
            Designed to work with matplotlib. It returns the starting coordinates
            and ending coordinates for a straight line that corresponds to a bond.
            When you call the object do this:
            x,y,z = object.get_bonds()
            axis3D.plot(x,y,z)
            This will give you a bond. Only works with water molecules created
            using the .water() method right now. Might extend it to simple 
            molecules with single atom centres.

        """
        return [([i[0][0],i[1][0]],[i[0][1],i[1][1]],[i[0][2],i[1][2]]) for i in self.bonds]
    def _2dfield(self, **kwargs):
        """
        This is not at all generalised. I just made this to make a 
        quick 2d vector plot.

        Returns
        -------
        Vector field parameters
        """  
        coords = (lambda a, b, c: ([[(a.append(i), b.append(j)) for j in c] for i in np.linspace(0, 3, 5)], (a, b)))([], [], np.linspace(.5, kwargs.get("bounds", 2.5), kwargs.get("res", 5)))[1]
        xout, zout = [], []
        for i in range(len(coords[0])):
            field = self.electric_field([coords[0][i], 0, coords[1][i]])
            vector = field[0] * field[1]
            xout.append(vector[0]), zout.append(vector[2])
        return coords[0], coords[1], xout, zout
            
        return coords
    def vector_field(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs : TYPE
            Includes: Size, resolution, centre
            The names are self explanatory

        Returns
        -------
        positions : list of arrays
            the positoin in space at which the electric field was determined
        vectors : list of tuples of arrays
            the electric field vector at that point in space 

        """
        size = kwargs.get("size", 3)
        resolution = kwargs.get("resolution", 5)
        scale = kwargs.get("vector_scale", 1e-10)
        centre = np.array(kwargs.get("centre", [0,0,0]))
        coords = np.linspace(-size, size, resolution + 1).tolist()
        #positions = []
        if kwargs.get("surface", False):
            positions = [np.array((i,j,4))+ centre for j in coords for i in coords]
        else:
            positions = [np.array((i,j,k))+ centre for k in coords for j in coords for i  in coords]
        vectors = [self.electric_field(pos) for pos in positions]
        return positions,  [i[0]*i[1]*scale for i in vectors]

def plot_charges(system, **kwargs):
    """
    This function is for plotting charges, bonds and vector fields in 3d space
    using matplotlib Axes3D
    Parameters
    ----------
    system : charge_sytem object
        charge_system object
    **kwargs : various types
            The name says what it does
            positions and vecors keywords pertain to the E field vectors and 
            their position in space, used to plot quivers.
    Returns
    -------
    3D plot of your system with the given parameters

    """
    pos, vec = kwargs.get("positions", []), kwargs.get("vectors", [])
    s = kwargs.get("s", 3)
    a = kwargs.get("a", s)
    el, az = kwargs.get("elevation", 30), kwargs.get("azimuth", -30)
    atom_size = kwargs.get("atom_size", 100)
    colour = kwargs.get("colour", "green")
    plt.figure()
    ax = plt.subplot(1,1,1, projection = "3d")
    scale = "%.2e" % system.position_factor
    ax.set_xlabel("x " + scale + " m"), ax.set_ylabel("y " + scale + " m"), ax.set_zlabel("z " + scale + " m")
    ax.set_title(kwargs.get("title", " "))
    ax.set_xlim(-s, s), ax.set_ylim(-s, s), ax.set_zlim(-s, a)  
    [ax.scatter3D(i.position[0], i.position[1], i.position[2], color = i.colour, s=atom_size) for i in system.charges] 
    for i in system.get_bonds():
        x,y,z = i
        ax.plot3D(x,y,z, color = "grey")
    for i in range(len(pos)):
        ax.quiver(pos[i][0],pos[i][1] ,pos[i][2], vec[i][0],vec[i][1],vec[i][2], color = colour)
    ax.view_init(elev = el, azim = az)

