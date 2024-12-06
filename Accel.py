# The acceleration event.
# Calculates the load on each wheel at each time step and
# ues motor torque and rpm to generate a longitudinal force
# from the tire model equation.

# Checklist of things to include
#   Motor efficiency curve
#   Pedal Map
#   Power limit
#   TCTV
#   Driver CG changes **
#   DRS
import math
import matplotlib.pyplot as plt
import numpy as np


class Accel:

    def __init__(self):
        # Track Length
        self.length = 75.02  # m
        #  Time Step
        self.delt = 0.001  # s
        # Car Speed
        self.vel = 0  # m/s
        # Car Acceleration
        self.accel = 0  # m/s^2
        # Current Load on Rear Tires
        self.rLoad = 0  # N
        # Current Load on Front Tires
        self.fLoad = 0  # N
        # Max usable power for motors
        self.usablePower = 72000  # W

        # For Graphs
        self.t = []  # time
        self.a = []  # acceleration
        self.v = []  # velocity
        self.d = []  # distance
        self.rl = []  # rear load
        self.fl = []  # front load
        self.rtq = []  # rear torque
        self.ftq = []  # front torque
        self.rpm = []  # rpm
        self.p = []  # power

    def run(self, car, graph, mode):

        self.t = []  # time
        self.a = []  # acceleration
        self.v = []  # velocity
        self.d = []  # distance
        self.rl = []  # rear load
        self.fl = []  # front load
        self.rtq = []  # rear torque
        self.ftq = []  # front torque
        self.rpm = []  # rpm
        self.p = []  # power

        self.vel = 0
        self.accel = 0
        vmax = car.redline * (2 * math.pi * car.RR) / (car.gearRatio * 60)

        time = 0
        dist = 0

        while dist < self.length:
            dist += self.vel * self.delt
            self.accelerate(car, vmax, mode)

            # For Graphs
            self.t.append(time)
            self.d.append(dist)
            self.v.append(self.vel)
            self.a.append(self.accel)
            self.rl.append(self.rLoad)
            self.fl.append(self.fLoad)

            time += self.delt

        print("Completed in", round(time, 3), "seconds")
        if graph:
            self.graph()

        return [self.d, self.v]

    def accelerate(self, car, vmax, mode):
        self.vel += self.accel * self.delt  # add velocity
        if self.vel > vmax:  # check velocity is not greater than top speed
            self.vel = vmax

        rpm = self.vel * 60 * car.gearRatio / (car.RR * 2 * math.pi)  # get current rpm
        maxT = car.maxT(rpm)

        self.cornerLoads(car)  # get loads on each corner in this instant
        # optimal power distribution mode (accounts motor efficiency)
        if mode == 1:
            # finding maximum possible torque to wheels
            maxT_f = car.getTorqueLimit(self.fLoad / 2)
            T_f = maxT_f if maxT_f < maxT else maxT
            maxT_r = car.getTorqueLimit(self.rLoad / 2)
            T_r = maxT_r if maxT_r < maxT else maxT

            if rpm < 100:
                T_f -= 2
            # check if max torque exceeds power limit
            powerLimited = True if 2 * (car.motorPower(T_f, rpm) + car.motorPower(T_r, rpm)) > self.usablePower else False

            if powerLimited:
                sol = self.distributePower(T_f, T_r, rpm, car)  # finds a power distribution that results in max total torque
                frontT = min(sol)
                rearT = max(sol)
            else:
                frontT = T_f
                rearT = T_r
        # lame power bias to rear way
        elif mode == 2:

            maxT_r = car.getTorqueLimit(self.rLoad / 2)
            if maxT_r > maxT:
                maxT_r = maxT
            maxPower_r = 2 * car.motorPower(maxT_r, rpm)
            if maxPower_r > self.usablePower:
                rearT = car.powerToTorque(self.usablePower, rpm)
                frontT = 0
            else:
                rearT = maxT_r
                maxT_f = car.getTorqueLimit(self.fLoad / 2)
                frontT = car.powerToTorque(self.usablePower - maxPower_r, rpm)
                if frontT > maxT_f:
                    frontT = maxT_f
                if frontT > maxT:
                    frontT = maxT
                if (rpm < 100):
                    frontT = 19
        #net force from torque - aero drag
        force = ((2 * car.gearRatio * (rearT + frontT)) / car.RR) - (car.cld * self.vel**2)
        self.accel = force / car.mass

        # For Graphs
        self.rtq.append(rearT)
        self.ftq.append(frontT)
        self.rpm.append(rpm)
        self.p.append(2 * (car.motorPower(rearT, rpm) + car.motorPower(frontT, rpm)))

    def cornerLoads(self, car):

        downForce = car.clv * self.vel**2

        # Longitudinal Non-Suspended Load Transfer
        NS_LT = (car.unsprungMass * self.accel * car.unsprungCG[2]) / car.wheelBase
        # Longitudinal Suspended Geometric Load Transfer
        SG_LT = (car.sprungMass * self.accel * car.PC[2]) / car.wheelBase
        # Longitudinal Suspended Elastic Load Transfer
        SE_LT = (car.sprungMass * self.accel * (car.sprungCG[2]-car.PC[2])) / car.wheelBase
        # Total Load Transfer
        loadTransfer = NS_LT + SG_LT + SE_LT

        self.fLoad = car.fLoad - loadTransfer + downForce * car.aerobias
        self.rLoad = car.rLoad + loadTransfer + downForce * (1-car.aerobias)

    # Finds the max torque distribution and a given rpm
    # ensures torque distribution will not slip the wheels
    def distributePower(self, T_f, T_r, rpm, car):
        aprx = (self.usablePower / 4) * 60 / (2 * math.pi * rpm)  # guess even distribution
        if rpm > 4000 and rpm < 18000:  # error adjustment for faster code
            error = 35000/(rpm-2700) - 2
            aprx -= error

        if aprx > T_f:
            aprx = T_f
        end = aprx - 5
        if end < 0:
            end = 0
        range = np.arange(aprx, end, -0.1)  # testing range
        max = 0
        sol = [0, 0]
        #print(T_f, T_r, rpm, self.d[-1])
        for i in range:
            tq = car.powerToTorque((self.usablePower - 2 * car.motorPower(i, rpm)), rpm)
            if tq > T_r:
                tq = T_r
            total = i + tq
            if total > max:
                max = total
                sol = [i, tq]
            elif max * 0.995 > total:
                return sol
        return sol

    def graph(self):
        fig, axs = plt.subplots(3, figsize=(12, 10))
        fig.suptitle('Motion Graphs')
        axs[0].plot(self.t, self.d)
        axs[0].set_title('Distance v Time')
        axs[1].plot(self.t, self.v)
        axs[1].set_title('Velocity v Time')
        axs[2].plot(self.t, self.a)
        axs[2].set_title('Acceleration v Time')

        plt.figure(2, figsize=(12, 10))
        plt.plot(self.t, self.rl, color='r', label='rear')
        plt.plot(self.t, self.fl, color='b', label='front')
        plt.xlabel('Time')
        plt.ylabel('Load (N)')
        plt.title('Load v Time')
        plt.legend()

        fig, axs = plt.subplots(3, figsize=(12, 10))
        fig.suptitle('Motor Graphs')
        axs[0].plot(self.t, self.rtq, color='r', label='rear')
        axs[0].plot(self.t, self.ftq, color='b', label='front')
        axs[0].set_title('Torque v Time')
        axs[0].legend()
        axs[1].plot(self.t, self.rpm)
        axs[1].set_title('RPM v Time')
        axs[2].plot(self.t, self.p)
        axs[2].set_title('Power v Time')

        plt.show()

