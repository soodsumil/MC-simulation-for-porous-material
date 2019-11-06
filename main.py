import numpy as np 
import penvs 
import physics 
import particles 
import random
import pylab 
import geom

## This defines a rectangular bed as the boundary of the porous media. xmin=-100, xmax=100, ymin=0, ymax=100. xmin and xmax are not very immportant, because the 
## structures is considered periodic, i.e if a ray passes x>100, it comes bac from the other end x=-100.

def example():
    B = penvs.bed_boundary(-100,100,0,100)

    els = [particles.ellipse_particle(x=B.random_x(margin=10), 
                                     y=B.random_y(margin=10),
                                     r1=random.uniform(1,10),
                                     r2=random.uniform(1,10),
                                     alpha=random.uniform(0,2*np.pi))
                                     for i in range(100)]
    E = penvs.distributed_ellipse(elipses=els, 
                                 medium_particle=penvs.medium(0.8,1.5e-8),
                                 medium_void=penvs.medium(1.0,0.0),
                                 bndry=B)
                            
    r = physics.ray(wlen=1.0e-6,
                    env=E,pos=np.array([0,0,]),
                    direction=geom.vec_normalize(np.array([random.uniform(-1,1),random.uniform(0.1,1.0)])),
                    polar=None
                    )
   
    ## This perforoms full cycle of ray tracing until the ray is outside the medium boundaries (B) or its power is smaller than the tolerance
    r.full_cycle()
    r.draw()
    pylab.show()

if __name__ == "__main__":
    example()



