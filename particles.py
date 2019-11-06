from numpy import array,sin,cos,linspace,pi,dot,sqrt,vstack,linalg
import pylab as plt 
from math import asin
import geom
class ellipse_particle(object): 
    '''
        defines a 2D ellipse partciel with center at (x,y), radii=(r1,r2) and rotation (alpha). 
    '''
    def __init__(self,x,y,r1,r2,alpha):
        self.center_ = array([x,y])
        self.radii_ = array([r1,r2])
        self.rotation_alpha_ = alpha 
        self.u_ = array([cos(alpha), sin(alpha)])
        self.v_ = array([-sin(alpha), cos(alpha)])
    def boundary_point(self,theta):
        return self.center_ + self.radii_[0]*cos(theta)*self.u_ + self.radii_[1]*sin(theta)*self.v_ 
    def draw(self): 
        x = array([ list( self.boundary_point(t) ) for t in linspace(0,2*pi,1000) ])
        plt.plot(x[:,0],x[:,1],'-') 
        #plt.show()
    def center_distance(self,x):
        return sqrt( dot( x-self.center_ , self.u_ )**2 / self.radii_[0]**2  +  dot( x-self.center_ , self.v_ )**2 / self.radii_[1]**2)

    def is_inside(self,x):
        return self.center_distance(x) <= 1

    def tangent_direction_on_boundary(self,theta,eps=0.01*pi):
        x1 = self.boundary_point(theta-eps)
        x2 = self.boundary_point(theta+eps)
        return x1 - x2

    def segment_intersect(self,p1,p2,tol=1.0e-6):
        p_inside = None
        p_outside = None
        if self.is_inside(p1): p_inside = p1 
        else: p_outside = p1 
        if self.is_inside(p2): p_inside = p2
        else: p_outside = p2 
        
        if (p_inside is None) or (p_outside is None): 
            raise AssertionError('one of p1 and p2 must be inside and the other must be outside the ellipse')
        while geom.vec_norm(p_inside-p_outside) >= tol*sum(self.radii_) : 
            p = 0.5*(p_inside + p_outside)
            if self.is_inside(p): p_inside = p 
            else: p_outside = p 

        return p_inside 
    
    def theta_of_bounary_point(self,p,tol=1.0e-3):
        assert abs(1.0 - self.center_distance(p)) <= tol*sum(self.radii_), 'point p not on boundary'
        H = vstack([self.radii_[0]*self.u_, self.radii_[1]*self.v_]).T
        cs = dot(linalg.inv(H),(p - self.center_))
        t = asin(cs[1])
        if cs[0] < 0:
            t = (pi - t)  
        return t
        

         
    
