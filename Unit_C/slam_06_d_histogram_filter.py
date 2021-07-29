# Histogram implementation of a bayes filter - combines
# convolution and multiplication of distributions, for the
# movement and measurement steps.
# 06_d_histogram_filter
# Claus Brenner, 28 NOV 2012
from pylab import plot, show, ylim
from distribution import *

def move(distribution, delta):
    """Returns a Distribution that has been moved (x-axis) by the amount of
       delta."""
    return Distribution(distribution.offset + delta, distribution.values)



# --->>> Copy your convolve(a, b) and multiply(a, b) functions here.

def multiply(a, b):
    """Multiply two distributions and return the resulting distribution."""

    # --->>> Put your code here.
    start_a = a.start()
    end_a = a.stop()
    start_b = b.start()
    end_b = b.stop()
    
    start = min(start_a, start_b)
    end = max(end_a, end_b)
    values = []
    for i in range(start, end):
        values.append(a.value(i) * b.value(i))
    
    a = Distribution(start, values)
    a.normalize()
        
    return a  # Modify this to return your result.

def convolve(a, b):
    """Convolve distribution a and b and return the resulting new distribution."""

    # store the previous offset
    new_offset = a.offset + b.offset
    # initialize an array for the posterior distribution of size: size(pos)+2
    new_values = [0] * ((len(b.values)+(len(a.values)-1)))
    # outter for loop, loops over the prior distribution (size: len(move))
    for i in range(len(a.values)):
    # inner for loop, loops over the individual conditional probability (size: len(prior))
        for j in range(len(b.values)):
            new_values[i+j] += a.values[i]*b.values[j]
        
    return Distribution(new_offset, new_values)  # Replace this by your own result.


if __name__ == '__main__':
    arena = (0,220)

    # Start position. Exactly known - a unit pulse.
    start_position = 10
    position = Distribution.unit_pulse(start_position)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         linestyle='-', drawstyle = 'steps')

    # Movement data.
    controls  =    [ 20 ] * 10

    # Measurement data. Assume (for now) that the measurement data
    # is correct. - This code just builds a cumulative list of the controls,
    # plus the start position.
    p = start_position
    measurements = []
    for c in controls:
        p += c
        measurements.append(p)

    # This is the filter loop.
    for i in range(len(controls)):
        # Move, by convolution. Also termed "prediction".
        control = Distribution.triangle(controls[i], 10)
        position = convolve(position, control)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='b', linestyle='-', drawstyle = 'steps')

        # Measure, by multiplication. Also termed "correction".
        measurement = Distribution.triangle(measurements[i], 10)
        position = multiply(position, measurement)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='r', linestyle='-', drawstyle = 'steps')

    show()
