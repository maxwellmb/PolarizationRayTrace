import numpy as np


def get_uncoated_amplitude_coefficients(aoi,n0,n1,reflected=True):
    """
    Return the reflection or transmission amplitude coefficients
    for an uncoated surface

    Args:
        aoi ([float32]): Angle of Incidence in radians
        n0 ([type]): Current index of refraction
        n1 ([type]): Index of refraction of the surface
        reflection (bool, optional): Toggle the output between "reflection" 
                            and "transmission". Defaults to True.
    """
    n = n1/n0

    if reflected:
        #The s-polarization - CLY Equation 8.12
        r_s = (np.cos(aoi) - np.sqrt(n**2-np.sin(aoi)**2))/ \
              (np.cos(aoi) + np.sqrt(n**2-np.sin(aoi)**2))
        
        #The p-polarization - CLY Equation 8.14
        r_p = (n**2*np.cos(aoi) - np.sqrt(n**2-np.sin(aoi)**2))/ \
              (n**2*np.cos(aoi) + np.sqrt(n**2-np.sin(aoi)**2))
        return r_s,r_p
    else: #Else - transmission
        #The s-polarization - CLY Equation 8.13
        t_s = 2*np.cos(aoi)/(np.cos(aoi)+np.sqrt(n**2-np.sin(aoi)**2))
        #The p-polarization - CLY Equation 8.15
        t_p = 2*n*np.cos(aoi)/(n**2*np.cos(aoi)+np.sqrt(n**2-np.sin(aoi)**2))
        return t_s,t_p

def get_single_film_amplitude_coefficients(aoi,n0,n1,n2,layer_thickness, 
                                            wavelength,reflected=True):
    """Return the reflection or transmission amplitude coefficients
    for a surface with a thin film coating. See Section 13.2 of CLY

    Args:
        aoi ([float32]): Angle of Incidence in radians
        n0 ([float32]): Current index of refraction (maybe complex)
        n1 ([float32]): Index of refraction of the thin film (maybe complex)
        n2 ([float32]): Index of refraction of the main surface (maybe complex)
        layer_thickness ([float32]): thickness of thin film in m
        wavelength ([float32]): wavelength of light in m
        reflection (bool, optional): Toggle the output between "reflection" 
                            and "transmission". Defaults to True.
    """    


    cos0 = np.cos(aoi)
    cos1 = np.sqrt(1-(n0**2*np.sin(aoi)**2)/n1**2)
    cos2 = np.sqrt(1-(n0**2*np.sin(aoi)**2)/n2**2)

    #CLY Equation 13.12
    phase_thickness = 2*np.pi/wavelength*n1*layer_thickness*cos1

    if reflected: 
        #CLY Equations 13.4-13.7
        r01_s = (n0*cos0-n1*cos1)/(n0*cos0+n1*cos1)
        r10_s = - r01_s
        r01_p = (n1*cos0-n0*cos1)/(n1*cos0+n0*cos1)
        r10_p = - r01_p

        r12_s = (n1*cos1-n2*cos2)/(n1*cos1+n2*cos2)
        r12_p = (n2*cos1-n1*cos2)/(n2*cos1+n1*cos2)

        #CLY Equation 13.13
        r_s = (r01_s+r12_s*np.exp(2j*phase_thickness))/(1-r12_s*r10_s*np.exp(2j*phase_thickness))
        #CLY Equation 13.14
        r_p = (r01_p+r12_p*np.exp(2j*phase_thickness))/(1-r12_p*r10_p*np.exp(2j*phase_thickness))

        return r_s,r_p
    else: 
        #CLY Equations 13.8-13.11
        t01_s = 2*n0*cos0/(n0*cos0+n1*cos1)
        t01_p = 2*n0*cos0/(n1*cos0+n0*cos1)

        t12_s = 2*n1*cos1/(n1*cos1+n2*cos2)
        t12_p = 2*n1*cos1/(n2*cos1+n1*cos2)

        t_s = t01_s*t12_s/(1-r01_s*r12_s*np.exp(2j*phase_thickness))
        t_p = t01_p*t12_p/(1-r01_p*t12_p*np.exp(2j*phase_thickness))

        return t_s,t_p

