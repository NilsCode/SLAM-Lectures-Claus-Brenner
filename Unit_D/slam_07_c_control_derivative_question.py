# This adds the derivative of g, this time with respect to the control
# (left and right motor movement).
#
# slam_07_c_control_derivative
# Claus Brenner, 11.12.2012
from lego_robot import *
from math import sin, cos, pi
from numpy import *

class ExtendedKalmanFilter:

    @staticmethod
    def g(state, control, w):
        x, y, theta = state
        l, r = control
        if r != l:
            alpha = (r - l) / w
            rad = l/alpha
            g1 = x + (rad + w/2.)*(sin(theta+alpha) - sin(theta))
            g2 = y + (rad + w/2.)*(-cos(theta+alpha) + cos(theta))
            g3 = (theta + alpha + pi) % (2*pi) - pi
        else:
            g1 = x + l * cos(theta)
            g2 = y + l * sin(theta)
            g3 = theta

        return array([g1, g2, g3])

    @staticmethod
    def dg_dcontrol(state, control, w):
        theta = state[2]
        l, r = tuple(control)
        alpha = (r - l) / w
        rad = l/alpha
        m = zeros([3,2])
        
        if r != l:
            theta_ = theta + alpha
            m[0,0] = w*r*(sin(theta_) - sin(theta))/((r - l)**2) - (r + l)*cos(theta_)/(2*(r-l))
            m[1,0] = w*r*(-cos(theta_) + cos(theta))/((r - l)**2) - (r + l)*sin(theta_)/(2*(r-l))
            m[2,0] = -1/w
            
            m[0,1] = -w*l*(sin(theta_) - sin(theta))/((r - l)**2) + (r + l)*cos(theta_)/(2*(r-l))
            m[1,1] = -w*l*(-cos(theta_) + cos(theta))/((r - l)**2) + (r + l)*sin(theta_)/(2*(r-l))
            m[2,1] = 1/w
            
        else:
            theta_ = theta
            m[0,0] = 0.5*(cos(theta) + l*sin(theta)/w)
            m[1,0] = 0.5*(sin(theta) - l*cos(theta)/w)
            m[2,0] = -1/w
            
            m[0,1] = 0.5*(-l*sin(theta)/w + cos(theta))
            m[1,1] = 0.5*(l*cos(theta)/2 + sin(theta))
            m[2,1] = 1/w
            
        return m


if __name__ == '__main__':
    # If the partial derivative with respect to l and r (the control)
    # are correct, then the numerical derivative and the analytical
    # derivative should be the same.

    # Set some variables. Try other variables as well.
    # In particular, you should check cases with l == r and l != r.
    x = 10.0
    y = 20.0
    theta = 35. / 180. * pi
    state = array([x, y, theta])
    l = 50.0
    r = 54.32
    control = array([l, r])
    w = 150.0

    # Compute derivative numerically.
    print( "Numeric differentiation dl, dr")
    delta = 1e-7
    control_l = array([l + delta, r])
    control_r = array([l, r + delta])
    dg_dl = (ExtendedKalmanFilter.g(state, control_l, w) -\
             ExtendedKalmanFilter.g(state, control, w)) / delta
    dg_dr = (ExtendedKalmanFilter.g(state, control_r, w) -\
             ExtendedKalmanFilter.g(state, control, w)) / delta
    dg_dcontrol_numeric = column_stack([dg_dl, dg_dr])
    print( dg_dcontrol_numeric)

    # Use the above code to compute the derivative analytically.
    print( "Analytic differentiation dl, dr:")
    dg_dcontrol_analytic = ExtendedKalmanFilter.dg_dcontrol(state, control, w)
    print( dg_dcontrol_analytic)

    # The difference should be close to zero (depending on the setting of
    # delta, above).
    print( "Difference:")
    print( dg_dcontrol_numeric - dg_dcontrol_analytic)
    print( "Seems correct:", allclose(dg_dcontrol_numeric, dg_dcontrol_analytic))
