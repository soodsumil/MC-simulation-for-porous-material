from numpy import array,sqrt,arctan,pi,cos,sin
from math import asin 
def vec_norm(x): 
    return sqrt( sum(x**2) )
def vec_normalize(x):
    return x/vec_norm(x) 

#def normal_axis(tang):
    ## counter clockwise normal to tang vector:
    #return vec_normalize(array([tang[1],tang[0]]))

def vec_orthog(x):
    return vec_rotate(x,0.5*pi)

def angle_from_sin_cos(s,c,tol=1.0e-3):
    assert abs(1.0-s*s-c*c) <= tol, 'sin^2 + cos^2 must be equal to 1.0'
    x = asin(s)
    if (c >= 0):
        return x 
    else:
        return pi - x


def vec_angle(x):
    return angle_from_sin_cos(x[1]/vec_norm(x), x[0]/vec_norm(x))

def angle_between(vec1,vec2):
    return vec_angle(vec1) - vec_angle(vec2)

def vec_rotate(x,angle):
    return array([cos(angle)*x[0] - sin(angle)*x[1], 
                  sin(angle)*x[0] + cos(angle)*x[1]
                  ])
#
def reflect_direction(incident,surface):
    theta = angle_between(surface,incident)
    return vec_rotate( incident, 2*theta )
#
#def angle_deviate(dir1, angle):
    
