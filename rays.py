import numpy as np
from interfaces import *
from surface import *

class Ray:
    """
    This class defines a ray

    Properties: 
    initial position
    propagation direction 

    current ray position
    current propagation direction

    wavelength
    """

    def __init__(self,r0,k0,n0 = None):
        

        self.r0 = r0 # The initial ray position
        self.k0 = k0 # The initial propagation direction

        if n0 is None:
            self.n0 = 1.000277 #The index of refraction
        else: 
            self.n0 = n0

        self.current_position = r0 #The current position
        self.current_propagation_direction =k0 #The current propagation direction
        self.current_n = self.n0 #The current index of refraction

        #And some lists to hold the histories of the ray
        self.past_positions = []
        self.past_directions = []
        self.past_n = [] 

        #The polarization matrices
        self.P_matrices = []

    def propagate_to_surface(self,surface,reflection=True):
        '''

        '''

        intercept,distance = surface.get_ray_intercept(self)

        self.past_positions.append(self.current_position)
        self.past_directions.append(self.current_propagation_direction)
        self.past_n.append(self.current_n)
        
        self.current_position = intercept
        if reflection:
            self.current_propagation_direction = surface.reflected_propagation_vector(self)
        else:
            self.current_propagation_direction = surface.refracted_propagation_vector(self)
    

        #If we're just encountering a plane, then we don't need to update the
        # index of refraction
        if not isinstance(surface,plane):
            # self.current_n = 
            
            if not reflection:
                self.current_n = surface.n_indices
        # else:
            # print("Gotcha!")


    def get_polarization_matrix(self,surface,wavelength=None,reflected=True):
        """[summary]

        Args:
            surface ([type]): [description]
            indices ([list of floats]): [description]
        """        

        n_indices = len(surface.n_indices)


        #The angle of incidence
        aoi = np.arccos(np.dot(self.past_directions[-1],surface.get_surface_normal(self.current_position)))
        
        #The AOI should be between 0 and 90
        if aoi > np.pi/2: 
            aoi = np.pi-aoi

        #Get the amplitude coefficients
        if n_indices == 1: 
            a_s, a_p = get_uncoated_amplitude_coefficients(aoi,
                            self.current_n,surface.n_indices[0],reflected=reflected)
        elif n_indices == 2:
            if wavelength is not None:
                a_s, a_p = get_single_film_amplitude_coefficients(aoi,
                        self.current_n,surface.n_indices[0],surface.n_indices[1],surface.layer_thicknesses[0],wavelength,reflected=reflected)
            else: 
                print("You need to pass in the layer thickess in meters, as a one element list to the keyword 'material_properties'")
                print("and a wavelength in meters to the 'wavelength' keyword")
                return
        else: 
            print("You've entered more than 2 indices, multilayer stacks of materials are not yet implemented")
            print("I'm not doing anything right now")
            return

        if reflected:
                local_polarization_matrix = np.array([[a_s,0,0],[0,-a_p,0],[0,0,1]])
        else:
                local_polarization_matrix = np.array([[a_s,0,0],[0,a_p,0],[0,0,1]])

        o_in,o_out = surface.calculate_orthogonal_matrices(self,reflected=reflected)

        # print(local_polarization_matrix)
        # print(o_in)

        # print(o_out)

        # return o_out@local_polarization_matrix@o_in

        #CLY Table 10.4
        #Python does the matrix multiplication in the opposite order here. 
        self.P_matrices.append(np.linalg.inv(o_in)@local_polarization_matrix@o_out)
        # self.P_matrices.append(o_out@local_polarization_matrix@np.linalg.inv(o_in))

    def get_jones_pupil(self, o_E,o_X):
        """[summary]

        Args:
            o_E (3x3 matrix): The orthogonal matrix to rotate the entrance 
                                pupil to the global coordinates
            o_X (3x3 matrix): The orthogonal matrix to rotate the exit pupil to the global coordinates
        """ 

        #First multiply all the P-matrices together
        
        #The number of P matrices after the first
        n_P = len(self.P_matrices)

        if n_P > 0:
            # P = self.P_matrices[-1]
            # if n_P > 1:
            #     for i in np.arange(2,n_P+1):
            #         #Multiply the matrices together
            #         P = P @  self.P_matrices[-i]

            P = self.P_matrices[0]
            if n_P > 1:
                for i in np.arange(1,n_P):
                    #Multiply the matrices together
                    P = P @  self.P_matrices[i]
            
            J = o_E @ P @ np.linalg.inv(o_X)

            #CLY - Equation 10.48
            # J = np.linalg.inv(o_X) @ P @ o_E

            self.J = J[:2,:2]


        else:
            self.J = np.identity(2)*np.nan

        
class Ray_Bundle:
    '''
    A bunch of Rays
    '''
    def __init__(self,rays = None):
        '''
        You can initiate without having any rays, or with a list or rays
        '''
        if rays is not None:
            self.rays = rays
        else: 
            self.rays = None
    
    def generate_ray_bundle(self,ranges, unit_vectors, plane_vertex, direction,nrays,n0=None,shape="square"):
        '''
        We'll generate a range of starting points in a plane, 
        but then project those on to the input unit_vectors for where we actually 
        want them in the global plane axes
        '''

        self.rays =[]

        if shape == "square":
            x_range = ranges[0]
            y_range = ranges[1]

            for i in np.linspace(x_range[0],x_range[1],nrays):
                for j in np.linspace(y_range[0],y_range[1],nrays):
                    location = i*np.array(unit_vectors[0]) + j*np.array(unit_vectors[1]) + plane_vertex
                    self.rays.append(Ray(location,direction,n0))
        
        # if shape == "polar":
        #     r_range = ranges[0]

        #     for 


    def propagate_to_surface(self,surface,reflection=True):
        '''
        Propagate all the rays to a given surface
        '''

        for ray in self.rays:
            if not np.isfinite(ray.current_position[0]):
                continue
                
            ray.propagate_to_surface(surface,reflection=reflection)

    def get_polarization_matrix(self,surface,wavelength=None,reflected=True):

        for ray in self.rays:
            if not np.isfinite(ray.current_position[0]):
                continue
            
            ray.get_polarization_matrix(surface,wavelength=wavelength,reflected=reflected)

    def get_jones_pupil(self,o_E):
        """
        Get all the jones pupils

        Args:
            o_E ([type]): [description]
        """        
        for ray in self.rays:

            #Double Dipole Coordinates - CLY Equation 11.11, 11.12
            k_exit = ray.current_propagation_direction
            x_exit = np.array([1-k_exit[0]**2/(1+k_exit[2]),-k_exit[0]*k_exit[1]/(1+k_exit[2]),-k_exit[0]])
            y_exit = np.array([-k_exit[0]*k_exit[1]/(1+k_exit[2]),1-k_exit[1]**2/(1+k_exit[2]),-k_exit[1]])

            #Let's instead put the anti-pole at a_loc=(0,1,0) - CLY Equation 11.14
            # x_exit = np.array([1-k_exit[0]**2/(1+k_exit[1]),-k_exit[0],-k_exit[0]*k_exit[2]/(1+k_exit[1])])
            # y_exit = np.array([k_exit[0]*k_exit[2]/(1+k_exit[1]),k_exit[2],k_exit[2]**2/(1+k_exit[1])-1])
            
            #Dipole Coordinates with poles at +/-x - CLT Equation 11.4
            # k_exit = ray.current_propagation_direction
            # y_exit = np.array([0,-k_exit[2],k_exit[1]])/np.sqrt(k_exit[2]**2+k_exit[1]**2)
            # x_exit = np.cross(y_exit,k_exit)
                
            o_X = np.array([x_exit,y_exit,k_exit])

            ray.get_jones_pupil(o_E,o_X)
        

