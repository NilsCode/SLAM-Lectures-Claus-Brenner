# The particle filter, prediciton and correction.
# In addition to the filtered particles, the density is estimated. In this
# simple case, the mean position and heading is computed from all particles.
#
# slam_08_c_density_estimation.
# Claus Brenner, 04.01.2013
from lego_robot import *
from slam_e_library import get_cylinders_from_scan, assign_cylinders
from math import sin, cos, pi, atan2, sqrt
import random
from scipy.stats import norm as normal_dist

import numpy as np

class ParticleFilter:
    def __init__(self, initial_particles,
                 robot_width, scanner_displacement,
                 control_motion_factor, control_turn_factor,
                 measurement_distance_stddev, measurement_angle_stddev):
        # The particles.
        self.particles = initial_particles

        # Some constants.
        self.robot_width = robot_width
        self.scanner_displacement = scanner_displacement
        self.control_motion_factor = control_motion_factor
        self.control_turn_factor = control_turn_factor
        self.measurement_distance_stddev = measurement_distance_stddev
        self.measurement_angle_stddev = measurement_angle_stddev

    # State transition. This is exactly the same method as in the Kalman filter.
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

        return (g1, g2, g3)

    def predict(self, control):
        """The prediction step of the particle filter."""
        left, right = control
        new_list = []
        # --->>> Put your code here.

        # Compute left and right variance.
        # alpha_1 is self.control_motion_factor.
        # alpha_2 is self.control_turn_factor.
        # Then, do a loop over all self.particles and construct a new
        # list of particles.
        # In the end, assign the new list of particles to self.particles.
        # For sampling, use random.gauss(mu, sigma). (Note sigma in this call
        # is the standard deviation, not the variance.)
        sig_left = (self.control_motion_factor*left)**2 + (self.control_turn_factor*(left-right))**2
        sig_right = (self.control_motion_factor*right)**2 + (self.control_turn_factor*(left-right))**2
        for i in range(len(self.particles)):
            l = random.gauss(left, sqrt(sig_left))
            r = random.gauss(right, sqrt(sig_right))
            
            new_state = self.g(self.particles[i], [l,r],self.robot_width)
            new_list.append(new_state)
        
        self.particles = new_list
        

    # Measurement. This is exactly the same method as in the Kalman filter.
    @staticmethod
    def h(state, landmark, scanner_displacement):
        """Takes a (x, y, theta) state and a (x, y) landmark, and returns the
           corresponding (range, bearing)."""
        dx = landmark[0] - (state[0] + scanner_displacement * cos(state[2]))
        dy = landmark[1] - (state[1] + scanner_displacement * sin(state[2]))
        r = sqrt(dx * dx + dy * dy)
        alpha = (atan2(dy, dx) - state[2] + pi) % (2*pi) - pi
        return (r, alpha)

    def probability_of_measurement(self, measurement, predicted_measurement):
        """Given a measurement and a predicted measurement, computes
           probability."""
        # Compute differences to real measurements.

        # --->>> Compute difference in distance and bearing angle.
        # Important: make sure the angle difference works correctly and does
        # not return values offset by 2 pi or - 2 pi.
        # You may use the following Gaussian PDF function:
        # scipy.stats.norm.pdf(x, mu, sigma). With the import in the header,
        # this is normal_dist.pdf(x, mu, sigma).
        # Note that the two parameters sigma_d and sigma_alpha discussed
        # in the lecture are self.measurement_distance_stddev and
        # self.measurement_angle_stddev.
        d_a = measurement[0]
        d_p = predicted_measurement[0]
        
        alpha_a = measurement[1]
        alpha_p = predicted_measurement[1]
        
        p_d = normal_dist.pdf(abs(d_a - d_p), 0, self.measurement_distance_stddev)
        angle = (alpha_a - alpha_p + pi) % (2*pi) - pi
        p_alpha = normal_dist.pdf(angle, 0, self.measurement_angle_stddev)
        
        p = p_d * p_alpha
        
        return p

    
    def compute_weights(self, cylinders, landmarks):
        """Computes one weight for each particle, returns list of weights."""
        weights = []
        for p in self.particles:
            # Get list of tuples:
            # [ ((range_0, bearing_0), (landmark_x, landmark_y)), ... ]
            assignment = assign_cylinders(cylinders, p,
                self.scanner_displacement, landmarks)

            # --->>> Insert code to compute weight for particle p here.
            # This will require a loop over all (measurement, landmark)
            # in assignment. Append weight to the list of weights.
            p_t = 1
            for tup in assignment:
                measurement = self.h(p, tup[1], self.scanner_displacement)
                predicted_measurement= tup[0]
                p_t *= self.probability_of_measurement(measurement, predicted_measurement)
                
            weights.append(p_t)  # Replace this.
        return weights

    def resample(self, weights):
        """Return a list of particles which have been resampled, proportional
           to the given weights."""

        # --->>> Insert your code here.
        # You may implement the 'resampling wheel' algorithm
        # described in the lecture.
        new_particles = []
        max_weigth = max(weights)
        M = len(weights)
        index = random.randint(0, len(weights) - 1)
        offset = 0
        for i in range(M):
            offset += random.uniform(0, 2*max_weigth)
            while offset >= weights[index]:
                index  = (index + 1) % M
                offset -= weights[index]
                
            new_particles.append(self.particles[index])
            
        return new_particles

    def correct(self, cylinders, landmarks):
        """The correction step of the particle filter."""
        # First compute all weights.
        weights = self.compute_weights(cylinders, landmarks)
        # Then resample, based on the weight array.
        self.particles = self.resample(weights)

    def print_particles(self, file_desc):
        """Prints particles to given file_desc output."""
        if not self.particles:
            return
        
        data_string = 'PA'
        #print >> file_desc, "PA",
        for p in self.particles:
            data_string += " %.0f %.0f %.3f" % p
            #print >> file_desc, "%.0f %.0f %.3f" % p,
        #print >> file_desc
        print(data_string, file = file_desc)
        
    def get_mean(self):
        """Compute mean position and heading from all particles."""

        # --->>> This is the new code you'll have to implement.
        # Return a tuple: (mean_x, mean_y, mean_heading).
        x_mean, y_mean, mh_x, mh_y,theta_mean = 0, 0 , 0, 0, 0
        for p in self.particles:
            x,y,theta = p
            x_mean += x
            y_mean += y
            mh_x += cos(theta)
            mh_y += sin(theta)
            
            
        x_mean /= len(self.particles)
        y_mean /= len(self.particles)
        mh_x /= len(self.particles)
        mh_y /= len(self.particles)
        theta_mean = atan2(mh_y, mh_x)
        return(x_mean , y_mean, theta_mean)
    '''
    def get_mean(self):
        """Compute mean position and heading from all particles."""

        # Return a tuple: (mean_x, mean_y, mean_heading).
        x = [ele[0] for ele in self.particles]
        y = [ele[1] for ele in self.particles]
        theta = [ele[2] for ele in self.particles]
        mean_x = sum(x)/len(x)
        meay_y = sum(y)/len(y)
        mean_heading_x = sum(np.cos(theta))
        mean_heading_y = sum(np.sin(theta))
        mean_heading = atan2(mean_heading_y, mean_heading_x)
        return (mean_x, meay_y, mean_heading)  # Replace this.
    '''


if __name__ == '__main__':
    # Robot constants.
    scanner_displacement = 30.0
    ticks_to_mm = 0.349
    robot_width = 155.0

    # Cylinder extraction and matching constants.
    minimum_valid_distance = 20.0
    depth_jump = 100.0
    cylinder_offset = 90.0

    # Filter constants.
    control_motion_factor = 0.35  # Error in motor control.
    control_turn_factor = 0.6  # Additional error due to slip when turning.
    measurement_distance_stddev = 200.0  # Distance measurement error of cylinders.
    measurement_angle_stddev = 15.0 / 180.0 * pi  # Angle measurement error.

    # Generate initial particles. Each particle is (x, y, theta).
    number_of_particles = 50
    measured_state = (1850.0, 1897.0, 213.0 / 180.0 * pi)
    standard_deviations = (100.0, 100.0, 10.0 / 180.0 * pi)
    initial_particles = []
    for i in range(number_of_particles):
        initial_particles.append(tuple([
            random.gauss(measured_state[j], standard_deviations[j])
            for j in range(3)]))

    # Setup filter.
    pf = ParticleFilter(initial_particles,
                        robot_width, scanner_displacement,
                        control_motion_factor, control_turn_factor,
                        measurement_distance_stddev,
                        measurement_angle_stddev)

    # Read data.
    logfile = LegoLogfile()
    logfile.read("robot4_motors.txt")
    logfile.read("robot4_scan.txt")
    logfile.read("robot_arena_landmarks.txt")
    reference_cylinders = [l[1:3] for l in logfile.landmarks]

    # Loop over all motor tick records.
    # This is the particle filter loop, with prediction and correction.
    f = open("particle_filter_mean.txt", "w")
    for i in range(len(logfile.motor_ticks)):
        # Prediction.
        control = map(lambda x: x * ticks_to_mm, logfile.motor_ticks[i])
        pf.predict(control)

        # Correction.
        cylinders = get_cylinders_from_scan(logfile.scan_data[i], depth_jump,
            minimum_valid_distance, cylinder_offset)
        pf.correct(cylinders, reference_cylinders)

        # Output particles.
        pf.print_particles(f)
        
        # Output state estimated from all particles.
        mean = pf.get_mean()
        print("F %.0f %.0f %.3f" %\
              (mean[0] + scanner_displacement * cos(mean[2]),
               mean[1] + scanner_displacement * sin(mean[2]),
               mean[2]), file = f)

    f.close()
