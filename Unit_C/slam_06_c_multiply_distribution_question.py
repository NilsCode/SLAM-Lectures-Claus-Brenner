# Multiply a distribution by another distribution.
# 06_c_multiply_distribution
# Claus Brenner, 26 NOV 2012
from pylab import plot, show
from distribution import *

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


if __name__ == '__main__':
    arena = (0,1000)

    # Here is our assumed position. Plotted in blue.
    position_value = 400
    position_error = 100
    position = Distribution.triangle(position_value, position_error)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         color='b', linestyle='-', drawstyle = 'steps')

    # Here is our measurement. Plotted in green.
    # That is what we read from the instrument.
    measured_value = 600
    measurement_error = 200
    measurement = Distribution.triangle(measured_value, measurement_error)
    plot(measurement.plotlists(*arena)[0], measurement.plotlists(*arena)[1],
         color='g', linestyle='-', drawstyle = 'steps')

    # Now, we integrate our sensor measurement. Result is plotted in red.
    position_after_measurement = multiply(position, measurement)
    plot(position_after_measurement.plotlists(*arena)[0],
         position_after_measurement.plotlists(*arena)[1],
         color='r', linestyle='-', drawstyle = 'steps')

    show()
