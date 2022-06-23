import rays
import surface

## Gemini Primary
gemini_primary_roc = 28.8 #meters
# Conic Constant 
gemini_primary_K = -1.003756
gemini_primary_diameter = 8.0\

X = [1,0,0]
Y = [0,1,0]
Z = [0,0,1]

gemini_primary = surface.conic_Surface([0,0,0],[Z,Y,X],gemini_primary_diameter,
                                       gemini_primary_roc,gemini_primary_K)


ray0_r0 = [0,0.,100]
ray0_k0 = [0,0,-1]

ray0 = rays.Ray(ray0_r0,ray0_k0)

ray0.propagate_to_surface(gemini_primary)

ray0.current_position