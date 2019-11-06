from numpy import cos,sin,sqrt,array
from math import asin
import geom
import pylab as plt


def reflection(theta_incident, n1, n2): 
    ti = theta_incident 
    if polar == 's': 
        return abs( (n1*cos(ti) - n2*sqrt(1-(n1*sin(ti)/n2)**2)) / 
                    (n1*cos(ti) + n2*sqrt(1-(n1*sin(ti)/n2)**2))
                    )**2
    elif polar  == 'p':
        return abs( (n1*sqrt(1-(n1*sin(ti)/n2)**2)-n2*cos(ti) ) / 
                    (n1*sqrt(1-(n1*sin(ti)/n2)**2)+n2*cos(ti) )
                    )**2
    else: 
        raise AssertionError("polar should be either 'p' or 's'")
    



def snell_angle(theta_incident, n1, n2):
    return asin(n1*sin(theta_incident)/n2)



def fresnel_reflection(theta,n1,n2,polar=None):
    assert cos(theta) >= 0
    if abs(n1/n2*sin(theta)) >= 1.0:
        Rs = 1
        Rp = 1
    else:
        Rs = abs((n1*cos(theta)-n2*sqrt(1-(n1/n2*sin(theta))**2)) / (n1*cos(theta)+n2*sqrt(1-(n1/n2*sin(theta))**2)) )
        Rp = abs((n1*sqrt(1-(n1/n2*sin(theta))**2)-(n2*cos(theta)))/(n1*sqrt(1-(n1/n2*sin(theta))**2)+(n2*cos(theta))))
    if polar == 'p':
        return Rp 
    elif polar == 's': 
        return Rs
    else: 
        return 0.5*(Rp + Rs)
    


class ray(object):
    ''' This defines a ray object. The ray should be associated an environment (env). It must be given an initial position (pos), which is a vector of length 2. It must be given a direction (direction), another vector of length 2, a polarization (polar) which is a string ('s' or 'p') or can be None for uniform polarization, an initial power (power), a step size (step) and a power tolerancce (power_tol). The ray is dead if its power is smaller than power_tol or is outside its environment.   
    '''
    def __init__(self,wlen,env,pos,direction,polar,power=1.0,step=0.1,power_tol=1.0e-5):
        self.pos_ = pos 
        self.direction_ = geom.vec_normalize( direction  )
        self.step_ = step 
        self.env = env 
        self.power_ = power 
        self.polar_ = polar 
        self.path_ = [pos]
        self.power_path_ = [power]
        self.power_tol_ = power_tol
        self.wlen_ = wlen
    def pos(self):
        return self.pos_
    def deltapos(self):
        return self.step_*self.direction_ 
    def record_(self):
        self.path_.append(self.pos_) 
        self.power_path_.append(self.power_)
    def move_(self):
        dpos = self.deltapos()
        self.power_ *= self.env.decay(self.pos_,dpos,self.wlen_)
        self.pos_ = self.pos_ + dpos
        self.pos_ = self.env.contain(self.pos_)
        self.record_()
    def draw(self,fancy=False):
        plt.subplot(2,1,1)
        if fancy:
            self.env.draw_fancy()
        else:
            self.env.draw()
        x = array(self.path_)
        #plt.plot(x[:,0],x[:,1],'k.',markersize=0.1)
        plt.plot(x[:,0],x[:,1],'k')
        plt.subplot(2,1,2)
        plt.plot(self.power_path_)
    
    def is_faded(self):
        return (self.power_ < self.power_tol_)
    def is_out(self):
        return self.env.is_out(self.pos_)
    def is_dead(self):
        return self.is_faded() or self.is_out()
    #def decay(self,beta,delta_s):
    #    self.power_ *= exp( -beta*delta_s)
    def traverse(self):
        if not self.is_dead():
            cpos = self.pos_ 
            dpos = self.deltapos()
            npos = cpos  + dpos  
            if self.env.media(cpos) == self.env.media(npos): 
                self.move_()
            else:
                #raise ValueError()
                self.direction_ = self.env.scatter(self.direction_,cpos,npos,self.polar_)
                self.move_()
    def full_cycle(self):
        while not(self.is_dead()):
            self.traverse()
