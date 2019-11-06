import physics 
import random 
import geom 
from numpy import *
import pylab as plt

class bed_boundary(object):
    def __init__(self,xmin,xmax,ymin,ymax,num_grid=100):
        self.xmin_ = xmin 
        self.xmax_ = xmax 
        self.ymin_ = ymin 
        self.ymax_ = ymax 
        self.num_grid_ = num_grid 
    def random_x(self,margin=0):
        return random.uniform(self.xmin_+margin,self.xmax_-margin)
    def random_y(self,margin=0):
        return random.uniform(self.ymin_+margin,self.ymax_-margin)
    def is_out(self,pos):
        return (pos[1] > self.ymax_) or (pos[1] < self.ymin_)
    def contain(self,pos):
        if (pos[0] > self.xmax_): return array([self.xmin_,pos[1]])
        if (pos[0] < self.xmin_): return array([self.xmax_,pos[1]])
        return pos 
    def draw(self):
        xspan = linspace(self.xmin_,self.xmax_,1000)
        yspan = linspace(self.ymin_,self.ymax_,1000)
        plt.plot(xspan,self.ymax_*ones(len(xspan)),'k-.')
        plt.plot(xspan,self.ymin_*ones(len(xspan)),'k-.')
        plt.plot(self.xmax_*ones(len(yspan)),yspan,'k-.')
        plt.plot(self.xmin_*ones(len(yspan)),yspan,'k-.')
    def iterate(self):
        for x in linspace(self.xmin_,self.xmax_,self.num_grid_):
            for y in linspace(self.ymin_,self.ymax_,self.num_grid_):
                yield array([x,y])


class medium(object):
    def __init__(self,n,k):
        self.k_ = k 
        self.n_ = n  
    def decay(self,dl,wlen):
        return exp(-2*pi*self.k_/wlen*dl )

class general_env(object):
    def media(self,pos):
        raise AssertionError('this method needs to be overridden')
        #return medium(0.0,1.0)
    def decay(self,pos,deltapos,wlen):
        return self.media(pos).decay(geom.vec_norm(deltapos),wlen)
    def auto_tangent(pos1,pos2):
        raise AssertionError('this method needs to be overridden')
    def is_out(self,pos):
        return self.bndry.is_out(pos)
    def contain(self,pos):
        return self.bndry.contain(pos)
    def scatter(self,direction,pos1,pos2,polar=None):
        tang  = self.auto_tangent(pos1,pos2)
        #ntang = geom.normal_axis(tang)  
        ntang = geom.vec_orthog(tang)
        #theta_incident = geom.angle_between(-direction, ntang) 
        theta_incident = geom.angle_between(ntang,-direction)
        if cos(theta_incident) < 0: 
            ntang = -ntang 
            theta_incident = geom.angle_between(ntang,-direction)

        n1 = self.media(pos1).n_ 
        n2 = self.media(pos2).n_

        R = physics.fresnel_reflection( theta_incident, n1, n2, polar=polar)
        u = random.uniform() 
        #raise ValueError()
        if u < R: 
            #return geom.vec_normalize( geom.reflect_direction(direction, ntang) )
            return geom.vec_normalize( geom.vec_rotate(-direction, 2*theta_incident) )
        else: 
            return geom.vec_normalize( geom.vec_rotate(-ntang, -physics.snell_angle(theta_incident,n1,n2) ) )
        #return array([1.0,0.0])
            

class distributed_ellipse(general_env):
    '''
       This class defines an environment which has rotated elipses in it. The inputs are a list of ellipses (elipses), medium_particle, medimum_void and boundary (bndry). medimum_particle and medimum_void must be instances of envs.medimum which take n,k values. n is the refractive index and k is for power dissipation. bndry is an instance of a bed_boundary object
    '''
    def __init__(self,elipses,medium_particle,medium_void, bndry):
        self.elipses_ = elipses 
        self.medium_particle_ = medium_particle 
        self.medium_void_ = medium_void 
        self.bndry =  bndry 
    def media(self,pos):
        inside_array = array([el.is_inside(pos) for el in self.elipses_])
        if any(inside_array): 
            return self.medium_particle_ 
        else: 
            return self.medium_void_ 
    def first_covering_ellipse(self,pos):
        inside_array = array([el.is_inside(pos) for el in self.elipses_])
        if any(inside_array): 
            i = where(inside_array)[0][0]
            return self.elipses_[i]
        else: 
            return None 

    def auto_tangent(self,pos1,pos2):
        el1 = self.first_covering_ellipse(pos1) 
        el2 = self.first_covering_ellipse(pos2)
        if not el1 is None:
            assert(el2 is None,'one of the points must be in void')
            el = el1 
        elif not el2 is None: 
            el = el2 
        else: 
            raise AssertionError('exactly one of the points must be covered by an ellipse')
        p = el.segment_intersect(pos1,pos2)
        t = el.theta_of_bounary_point(p)
        return el.tangent_direction_on_boundary(t)
                    

    def draw(self):
        self.bndry.draw()
        for el in self.elipses_: 
            el.draw()
    
    def draw_fancy(self):
        self.bndry.draw() 
        xs = []
        for x in self.bndry.iterate():
            if self.media(x) == self.medium_particle_: 
                xs.append(x)
        xs = array(xs)
        plt.plot(xs[:,0],xs[:,1],'r.')
            
        

    def porosity(self):
        s = 0; n = 0 
        for x in self.bndry.iterate():
            n += 1
            if self.media(x) == self.medium_void_: 
                s += 1 
        return float(s)/n
                

            
