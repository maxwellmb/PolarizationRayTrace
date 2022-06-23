## Most (if not all) of this is from the book: 
## Polarized Light and Optical Systems 
## by Russell A. Chipman, Wai-Sze Tiffany Lam and Garam Young
## Below I'll refer to this reference as CLY

import numpy as np
import copy
import scipy.optimize as so


class Surface:
    """
    A class used to represent a Surface

    based on Table 10.1 of Polarized Light and Optical Systems

    ...

    Attributes
    ----------
    vertex_location : [float,float,float]
        a 3D vector marking the location of the surface vectex
    surface_normal : [float, float, float]
        a 3D vector representing the surface normal at the vertex
    sound : str
        the sound that the animal makes
    num_legs : int
        the number of legs the animal has (default 4)

    Methods
    -------
    tbd(sound=None)
        Prints the animals name and what sound it makes

    aperture_equation(x,y,z)
        A binary map of the aperture extent

    """

    def __init__(self, vertex_location, surface_axes, diameter,n_indices = None,
                layer_thicknesses = None):
        '''
        Initialize
        '''
        #TODO: Some checks to make sure the inputs are the right format/dtype
        #
        self.vertex_location = vertex_location
        self.surface_axes = surface_axes
        self.diameter = diameter

        #If no index of refraction is provided, then assume air
        if n_indices is None: 
            self.n_indices = [1.000277]
        else: 
            #Make sure it's a list, even if it's only one element
            if not isinstance(n_indices, list):
                self.n_indices = [n_indices]
            else:
                self.n_indices = n_indices

        #Make sure it's a list, even if it's only one element
        if layer_thicknesses is not None:
            if not isinstance(layer_thicknesses, list):
                self.layer_thicknesses = [layer_thicknesses]
            else:
                self.layer_thicknesses = layer_thicknesses
        else: 
            self.layer_thicknesses = [10000] #If not provided make it really big


        ## The definition of the local axes will be an x,y,z axes and we will define the surface
        ## such that +ve x is the optical axis. 
        self.local_axes = np.array([[1,0,0],[0,1,0],[0,0,1]])

        ## We're going to later want to rotate things into the local frame we so
        ## we can just generate that rotation matrix right now. 
        ## We do this by contructing a matrix of dot-products of the surface normal in the 
        ## global coordinate axes with the unit vector in the local axes
        #         
        ## Take the outer product to get the rotation matrix: 

        self.rotation_matrix = np.tensordot(np.array(self.surface_axes),self.local_axes,axes=1)
    
    
    def surface_function(self,position):
        '''

        '''

        return None

    def aperture_function(self,position):
        '''
        Defines the 3D aperture for this surface 
        '''

        ### We'll define the aperture centered at the origin perpendicular 
        ### to the z-direction we'll rotate and shift the ray later. 


        x,y,_ = position
        
        if np.sqrt(x**2+y**2) <= self.diameter/2:
            return 1
        else:
            return 0

    def get_surface_normal(self,position):
        '''
        Get the surface normal at a given position
        '''

        return None

    def reflected_propagation_vector(self,ray):
        '''
        Get the reflected propagation vector for a Ray and a given position
        '''

        #The surface normal at the ray location
        surface_normal = self.get_surface_normal(ray.current_position)

        # print(surface_normal)
        #The reflected propagation direction in global coordinates- CLY Equation 10.17 
        k_new = ray.current_propagation_direction \
                - 2*(np.dot(ray.current_propagation_direction,surface_normal))*np.array(surface_normal)

        return k_new

    def refracted_propagation_vector(self,ray):
        '''
        Get the refracted propagation vector for a Ray and a given position
        '''

        #The surface normal at the ray location
        surface_normal = self.get_surface_normal(ray.current_position)

        #The incident angle on the surface - CLY Equation 10.16
        incident_angle = np.arccos(np.dot(ray.current_propagation_direction,surface_normal))

        #Cosine of the refracted angle - CLY Equation 10.19
        cos_refracted_angle = np.sqrt(1-(ray.current_n/self.n_indices[0])**2*(1-np.cos(incident_angle)**2))

        #The refracted propagation direction in global coordinates - CLY Equation 10.20
        k_new = ray.current_n/self.n_indices[0]*ray.current_propagation_direction \
                -(ray.current_n/self.n_indices[0]*np.cos(incident_angle)-cos_refracted_angle)*surface_normal
        
        return k_new

    def calculate_orthogonal_matrices(self,ray,reflected=True):
        '''
        Calculate the orthogonal matrices to rotate the incident basis into the exiting basis
        '''
        
        #The surface normal at the ray location
        surface_normal = self.get_surface_normal(ray.current_position)

        ### First calculate the incident s and p basis vectors: 
        #CLY Equation 10.22
        s_in = np.cross(ray.past_directions[-1],surface_normal)
        s_in /= np.linalg.norm(s_in)

        p_in = np.cross(ray.past_directions[-1],s_in)

        ### Calculate the incident orthogonal matrix (CLY - Table 10.4): 
        o_in = np.array([s_in,p_in,ray.past_directions[-1]])

        # print(surface_normal)
        # print(ray.current_propagation_direction)

        ### Calculate the outgoing s and p basis vectors: 
        #CLY Equation 10.22
        s_out = s_in

        #We calculate the new propataion vector when we propagate_to_surface
        k_out = ray.current_propagation_direction
        
        p_out = np.cross(k_out,s_out)

        ### Calculate the incident orthogonal matrix (CLY - Table 10.4): 
        if reflected:
            o_out = np.array([s_out,-p_out,k_out])
        else:
            o_out = np.array([s_out,p_out,k_out])
            

        return o_in,o_out

    def calculate_amplitude_coefficients(self,ray,position):
        '''
        '''
    
    def calculate_jones_matrix(self,ray,position):
        '''
        '''
    
    def calculate_polarization_matrix(self,ray,local_jones_matrix,reflected=True):
        '''
        Calculate the polarization ray tracing matrix - CLY Equation 10.31
        '''
        o_in, o_out = self.calculate_orthogonal_matrices(ray,reflected=reflected)

        j_extended = np.identity(3)
        j_extended[:2,:2] = local_jones_matrix 

        
    
    def calculate_parallel_transport_matrix(self,Ray,position):
        '''
        '''


class conic_Surface(Surface):
    '''
    '''
    def __init__(self,vertex_location,surface_normal,diameter,radius_of_curvature,
                conic_constant,side="top",n_indices=None,
                layer_thicknesses=None):
        super().__init__(vertex_location,surface_normal,diameter,n_indices=n_indices,
                            layer_thicknesses=layer_thicknesses)

        self.roc = radius_of_curvature
        self.K = conic_constant
        self.side = side # You could come at a conic section from the top of the bottom. 
                         # Top is the concave side 

    def get_ray_intercept(self,ray):
        '''
        Get the intercept of the input Ray for a conic surface
        '''

        ##### We want to find a position on the surface where the ray intercepts. 
        ### To do this we'll put the ray coordinates into the local coordinates of the 
        ### conic surface

        ## First we want to shift the zero-point from global to local coordinates
        ## (i.e. subtract the surface location in global coordinates from the ray position)
        ray_position = copy.deepcopy(ray.current_position)
        local_ray_position = np.array(ray_position) - np.array(self.vertex_location)

        ## We can apply a rotation matrix to rotate from the global coordites to the local
        ## Rotate the ray:
        # local_ray_position = np.matmul(self.rotation_matrix,local_ray_position)
        local_ray_position = self.rotation_matrix@local_ray_position
        
        ## Similarly we want to rotate the ray propagation vector
        
        local_ray_propagation_vector = self.rotation_matrix@ray.current_propagation_direction

        ## Here's a paper with instructions how to do it in 2D. Pretty simple. 
        # https://iopscience-iop-org.clsproxy.library.caltech.edu/article/10.1088/0959-5309/55/5/302/pdf

        # #It's not bad to extend this to 3D: it's still just a quadratic formula
        x0,y0,z0 = local_ray_position
        xv,yv,zv = local_ray_propagation_vector

        a = yv**2+zv**2+(self.K+1)*xv**2
        b = 2*y0*yv + 2*z0*zv - 2*self.roc*xv+2*(self.K+1)*x0*xv
        c = y0**2+z0**2-2*self.roc*x0+(self.K+1)*x0**2

        t1 = (-b + np.sqrt(b**2-4*a*c))/(2*a)
        t2 = (-b - np.sqrt(b**2-4*a*c))/(2*a)

        ## Because of how we've set up the conic, we want the one that will bring us closeset to the x-axis
        x_t1 = x0 + t1*xv
        x_t2 = x0 + t2*xv

        if np.abs(x_t1) < np.abs(x_t2):
            t = t1
        else: 
            t = t2

        ### Before I worked out the math I did it as below. Turns out that was more work that doing the math....
        # def to_minimize(t):

        #     #The x,y,z local coordinates after propagating a distance t
        #     x = local_ray_position[0] + t * local_ray_propagation_vector[0]
        #     y = local_ray_position[1] + t * local_ray_propagation_vector[1]
        #     z = local_ray_position[2] + t * local_ray_propagation_vector[2]

        #     # This is the radial coordinate often expressed just as why when writing
        #     # out the conic section formula in 2D with x as the optical axis. 
        #     yy = np.sqrt(y**2+z**2) 

        #     # Copying this formula from wikipedia: 
        #     # https://en.wikipedia.org/wiki/Conic_constant
        #     # but seems consistent with other sources
        #     conic_formula = yy**2 - 2*self.roc*x + (self.K+1)*x**2

        #     # if x < 0: #I don't want the "negative" hyperbola. 
        #     #     return 100000

        #     return np.abs(conic_formula)

        # #The initial guess will just be the current x-coordinate
        # t0 = np.abs(local_ray_position[0])
        # results = so.minimize(to_minimize,t0)
        # t = results.x[0]

        #The x,y,z local coordinates after propagating a distance t
        # x = local_ray_position[0] + t * local_ray_propagation_vector[0]
        # y = local_ray_position[1] + t * local_ray_propagation_vector[1]
        # z = local_ray_position[2] + t * local_ray_propagation_vector[2]

        # # This is the radial coordinate often expressed just as why when writing
        # # out the conic section formula in 2D with x as the optical axis. 
        # yy = np.sqrt(y**2+z**2) 

        #Now we can calculate the ray intercept position: 
        x_ray_intercept = local_ray_position[0] + t * local_ray_propagation_vector[0]
        y_ray_intercept = local_ray_position[1] + t * local_ray_propagation_vector[1]
        z_ray_intercept = local_ray_position[2] + t * local_ray_propagation_vector[2]

        ray_intercept = np.array([x_ray_intercept,y_ray_intercept,z_ray_intercept])

        #Check to see if we fall outside of the outer diameter
        radius = np.sqrt(ray_intercept[1]**2 + ray_intercept[2]**2)
        if radius > self.diameter/2:
            return [np.nan,np.nan,np.nan],0.

        ### Rotate back to global coordinates
        ray_intercept_global = np.linalg.solve(self.rotation_matrix,ray_intercept).T
        ray_intercept_global += np.array(self.vertex_location)

        return ray_intercept_global, t

    def get_surface_normal(self,position):
        '''
        Get the surface normal at a given position. 
        '''

        ### The input position should be in global coordinates. 
        ### We will first translate it to 

        position = copy.deepcopy(position)
        local_position = np.array(position) - np.array(self.vertex_location)

        ## We can apply a rotation matrix to rotate from the global coordites to the local
        ## Rotate the ray:
        # local_ray_position = np.matmul(self.rotation_matrix,local_ray_position)
        local_position = self.rotation_matrix@local_position
        
        ## To get the normal at this location we take the partial derivative of the 
        ## conic section equation for each variable: 

        normal = [-2*self.roc + 2*(self.K+1)*local_position[0],
                  2*local_position[1],
                  2*local_position[2]]

        #If we're using the botton side of the mirror, then flip the normal
        if self.side == "bottom":
            normal = -np.array(normal)

        
        #Normalize:
        normal /= np.linalg.norm(normal)
        
        #Convert back to global coordinates
        global_normal = np.linalg.solve(self.rotation_matrix,normal).T #rotate
        # global_normal += np.array(self.vertex_location) #Shift back to vertex location

        return global_normal

class plane(Surface):

    def __init__(self,vertex_location,surface_axes,diameter,surface_normal):
        super().__init__(vertex_location,surface_axes,diameter)

        self.surface_normal = surface_normal

    def get_ray_intercept(self,ray):
        '''
        Get the intercept of the input Ray with a plane
        '''
        #This wikipedia page gives the nice math that we want:
        # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection

        d = np.dot((self.vertex_location-ray.current_position), self.surface_normal)/ \
            np.dot(ray.current_propagation_direction, self.surface_normal)

        return ray.current_position + d*np.array(ray.current_propagation_direction),d

    def get_surface_normal(self,position):

        return self.surface_normal

    def reflected_propagation_vector(self,ray):
        return -ray.current_propagation_direction

    def refracted_propagation_vector(self,ray):
        return ray.current_propagation_direction

class flat_mirror(Surface):

    def __init__(self,vertex_location,surface_axes,diameter,surface_normal):
        super().__init__(vertex_location,surface_axes,diameter)

        self.surface_normal = surface_normal

    def get_ray_intercept(self,ray):
        '''
        Get the intercept of the input Ray with a plane
        '''
        #This wikipedia page gives the nice math that we want:
        # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection

        d = np.dot((self.vertex_location-ray.current_position), self.surface_normal)/ \
            np.dot(ray.current_propagation_direction, self.surface_normal)

        return ray.current_position + d*np.array(ray.current_propagation_direction),d

    def get_surface_normal(self,position):

        return self.surface_normal
