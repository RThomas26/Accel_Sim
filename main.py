# ACCEL SIM
# by Ruben Thomas 425-394-2271

# HOW TO USE
# The code runs an instance of the Car class in an instance of te Accel Class
# To run, first make an accel an excel object then call accel.run(Car, True/False, 1/2)
# specific details such as time step size and track length can be modified from within the accel class
# most commonly change time step. Smaller step is more accurate but also runs slower.
# No progress bar currently but if you go to the run method in accel and uncomment the print(dist) line
# you can gauge how far along each run is

# PARAMETERS
# Car()
# An instance of the Car class currently set to aprox T35
# Variables can be changed by either directly changing the car class
# or adding a parameter into the class construction and setting the desired
# variable to modify to that parameter. Useful for testing multiple setups in succession
# Currently below can be used as an example of how to run several accels testing different gear ratios
# True/False
# Boolean to chose whether graphs for  that specific accel run should be displayed after that run
# 1/2
# Chooses which power distribution method to use during accel
# 1 tries to find the optimal distribution such that the maximum amount of total torque is applied to the
# wheels without exceeding the power limit or breaking traction
# 2 simply puts the max amount of torque to the rears that will allow it to not slip and uses the remaining power
# to run the front wheels
# (option one is theoretically better but the graphs look funky so idk if it is bug free yet)

# Graphs
# If the True parameter is set graphs giving information about the torque, rpm, power, load dist, dist, vel, and accel
# all in regards to time are displayed
# No integrated and simplified way to do more specific graphs other than coding it yourself
# Get good at matplotlib
# An example of how I graphed the the distance v times of multiple
# accel runs on the same graph for comparison is below for reference

from Car import Car
from Accel import Accel
import time
import matplotlib.pyplot as plt
import numpy as np

start_time = time.time()

accel = Accel()
ratios = np.arange(12.5, 13.5, 0.2)


x = np.arange(0, 4, 0.001)

plt.figure(1, figsize=(12, 10))
plt.figure(2, figsize=(12, 10))

for r in ratios:
   print(r)
   result = accel.run(Car(r), True, 1)
   diff = len(x) - len(result[0])
   d = result[0] + ([0] * diff)
   v = result[1] + ([0] * diff)
   plt.figure(1)
   plt.plot(x, d, label=r)
   plt.figure(2)
   plt.plot(x, v, label=r)

plt.figure(1)
plt.xlabel('Time')
plt.ylabel('Dist')
plt.title('Dist v Time')
plt.legend()

plt.figure(2)
plt.xlabel('Time')
plt.ylabel('Vel')
plt.title('Vel v Time')
plt.legend()

plt.show()

end_time = time.time()
print("Run Time:", end_time-start_time, "seconds")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
