from __future__ import print_function
from astropy import units as u
from astropy import constants

"""
Examine the formation scenario of W51 from an energetics perspective

1. Assume there are two clouds, one at 55 km/s and one at 68 km/s.
They are 100's of pc apart, but are both interacting with a spiral arm.
The 68 km/s cloud has been accelerated from the tangent velocity of 54 km/s.

In this scenario, the 55 km/s cloud probably underwent a similar acceleration
in the past, but either self-collided or collided with other material and
decelerated to its current velocity.  The 68 km/s and 55 km/s clouds are not 
directly interacting.

Pros/Cons:
    + Explains the velocity of the 68 km/s cloud
    - Does not explain why the 55 km/s cloud is apprently younger than the
      W51 B 68 km/s cloud
    - Requires a deceleration mechanism for the 55 km/s cloud
"""

# 3.3e5 from co_intmaps 62-75 km/s
E68kms = (3.3e5 * u.M_sun * ((68-54)*u.km/u.s)**2).to(u.erg)
print("68 km/s cloud required {0} to reach current velocity.".format(E68kms))

"""
2. The W51 protocluster complex is the interaction point between the 50 km/s
cloud and the 68 km/s cloud.

 + Naturally explains the centroid velocity of the protocluster
 + Geometrically consistent along the line of sight
 0 Requires that the 50 km/s cloud at higher ell is moving toward the 68 km/s
   at lower ell.  If this motion is radial, it requires that both clouds are
   on the far side of the tangent point.
 - W51 B / 68 km/s appears older: if a spiral arm is responsible for triggering
   its star formation, it is moving in the wrong direction

"""

# Gravitational potential of the present-day protocluster
r = 15*u.pc
v = 13*u.km/u.s
m = (v**2 * r / constants.G).to(u.M_sun)
print("Enclosed mass required for a velocity {0} at distance "
      "{1}: {2}".format(r,v,m))

"""
3. The W51 cloud complex is purely hierarchical and turbulent, therefore there
is no instantaneous trigger

 + Natural the whole way through
 - Cannot explain the 68 km/s cloud's velocity
"""

"""
4. The W51 Main protocluster came from the 55 km/s cloud and has now expanded
in velocity by +/- 15 km/s, so any apparent interaction is illusory.

 - Cannot explain the 68 km/s cloud's velocity
"""

# again, enclosed mass here is key:
r = 5*u.pc
m = 2e5*u.M_sun
v = ((2*constants.G*m/r)**0.5).to(u.km/u.s)
print("delta-V = {0} for an enclosed mass of {1} at {2}".format(v,m,r))

"""
5. The W51 Main protocluster came from the 68 km/s cloud and has now expanded
in velocity to +/- 15 km/s.  It has decelerated to ~55-60 km/s.

 - Cannot explain the 68 km/s cloud's velocity
"""
