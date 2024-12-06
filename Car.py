# Car class contains all the relevant variables of the car in one place
import math

import sympy as sp
from sympy import symbols
from sympy.solvers import solve

from scipy.interpolate import griddata
from scipy.interpolate import CubicSpline

import pandas as pd

import numpy as np

import time


class Car:

    def __init__(self, var):
        self.carMass = 190.51 # kg (420 lb)
        self.sprungMass_NoDriver = 142.428 # kg (314 lb)
        self.unsprungMass = self.carMass - self.sprungMass_NoDriver
        self.driverMass = 54  # kg
        self.sprungMass = self.sprungMass_NoDriver + self.driverMass
        self.mass = self.unsprungMass + self.sprungMass

        self.sprungCG_NoDriver = [0.885, 0, 0.239]  # m
        self.unsprungCG = [0.299, 0, 0.087]  # m
        self.driverCG = [0.717, 0, 0.209]  # m
        self.sprungCG = self.combineCG(self.sprungCG_NoDriver, self.sprungMass_NoDriver, self.driverCG, self.driverMass)
        self.CG = self.combineCG(self.sprungCG, self.sprungMass, self.unsprungCG, self.unsprungMass)
        print(self.CG)

        self.frontalArea = 1.3657  # m^2
        # CHECK
        self.ClCd = 2.3
        self.CDA = 1.33
        self.Cl = 1
        #CHECK
        # Combo of all the constants before the velocity term for calcs
        self.clv = 3.287
        self.cld = 1.027
        self.aerobias = 0.71  # front bias

        # Suspension Coordinates
        # Taken from K-sketch
        # Origin is located in the middle of the front tire patches
        # Positive x is towards the back of the car
        # Positive y is towards right wheel
        # Positive z is up
        # All coordinates in meters

        # FR Suspension
        self.FR_Tire_Patch = [0, 0.624, 0]

        self.FR_Lower_Upright_Pickup = [-0.00462, 0.5988, 0.122]
        self.FR_Lower_Forward_Pickup = [-0.03, 0.21, 0.13]
        self.FR_Lower_Rear_Pickup = [0.12, 0.21, 0.13]

        self.FR_Upper_Upright_Pickup = [0.00586, 0.583, 0.298]
        self.FR_Upper_Forward_Pickup = [-0.10979, 0.232, 0.27634]
        self.FR_Upper_Rear_Pickup = [0.09043, 0.232, 0.26174]

        # RR Suspension
        self.RR_Tire_Patch = [1.545, 0.624, 0]

        self.RR_Lower_Upright_Pickup = [1.575, 0.595, 0.11]
        self.RR_Lower_Forward_Pickup = [1.49088, 0.27609, 0.1168]
        self.RR_Lower_Rear_Pickup = [1.62208, 0.27101, 0.11672]

        self.RR_Upper_Upright_Pickup = [1.545, 0.595, 0.305]
        self.RR_Upper_Forward_Pickup = [1.39508, 0.282, 0.2596]
        self.RR_Upper_Rear_Pickup = [1.63809, 0.27228, 0.272]

        # FL Suspension
        self.FL_Tire_Patch = [0, -0.624, 0]

        self.FL_Lower_Upright_Pickup = [-0.00462, -0.5988, 0.122]
        self.FL_Lower_Forward_Pickup = [-0.03, -0.21, 0.13]
        self.FL_Lower_Rear_Pickup = [0.12, -0.21, 0.13]

        self.FL_Upper_Upright_Pickup = [0.00586, -0.583, 0.298]
        self.FL_Upper_Forward_Pickup = [-0.10979, -0.232, 0.27634]
        self.FL_Upper_Rear_Pickup = [0.09043, -0.232, 0.26174]

        # RL Suspension
        self.RL_Tire_Patch = [1.545, -0.624, 0]

        self.RL_Lower_Upright_Pickup = [1.575, -0.595, 0.11]
        self.RL_Lower_Forward_Pickup = [1.49088, -0.27609, 0.1168]
        self.RL_Lower_Rear_Pickup = [1.62208, -0.27101, 0.11672]

        self.RL_Upper_Upright_Pickup = [1.545, -0.595, 0.305]
        self.RL_Upper_Forward_Pickup = [1.39508, -0.282, 0.2596]
        self.RL_Upper_Rear_Pickup = [1.63809, -0.27228, 0.272]

        self.wheelBase = self.RR_Tire_Patch[0] - self.FR_Tire_Patch[0]  # m
        # Weight Bias of Rear
        self.wbias = self.CG[0] / self.wheelBase
        # Load on a Rear Corner at Rest
        self.rLoad = self.mass * self.wbias * 9.81  # N
        # Load on a Front Corner at Rest
        self.fLoad = self.mass * (1 - self.wbias) * 9.81  # N

        # Longitudinal Only
        # Front IC
        self.frontIC = self.intersectFromPoints(self.RL_Lower_Forward_Pickup, self.RL_Lower_Rear_Pickup,
                                                self.RL_Lower_Upright_Pickup, self.RL_Upper_Forward_Pickup,
                                                self.RL_Upper_Rear_Pickup, self.RL_Upper_Upright_Pickup)
        # Rear IC
        self.rearIC = self.intersectFromPoints(self.FL_Lower_Forward_Pickup, self.FL_Lower_Rear_Pickup,
                                               self.FL_Lower_Upright_Pickup, self.FL_Upper_Forward_Pickup,
                                               self.FL_Upper_Rear_Pickup, self.FL_Upper_Upright_Pickup)
        # Pitch Center
        self.PC = self.intersectFromPoints(self.frontIC, self.RL_Tire_Patch, self.RL_Tire_Patch,
                                           self.rearIC, self.FL_Tire_Patch, self.FL_Tire_Patch)

        # Motor Efficiency
        self.motorPoints = None
        self.motorEff = self.createMotorData()
        self.gearRatio = var
        self.maxTorque = 21
        self.redline = 20000

        # Tire Model
        # Interpolated Eq of vertical load vs max longitudinal force
        self.tractionLimit = self.createTireTraction()  # N
        # Rolling Radius
        self.RR = 0.210  # m
        self.wheelDia = 0.4064  # m (16 in)

    # Returns the max longitudinal force based on the load on a tire
    def getTorqueLimit(self, load):
        return self.tractionLimit(load) * self.RR / self.gearRatio

    def maxT(self, rpm):
        lossrpm = 14620
        if (rpm > lossrpm):
            maxT = self.maxTorque - ((rpm - lossrpm) * ((self.maxTorque-15)/(self.redline-lossrpm)))
        else:
            maxT = self.maxTorque
        return maxT

    # Creates interpolated equation of load on tire to max longitudinal force
    def createTireTraction(self):
        df = pd.read_excel('Longitudinal_Force_Map.xlsx')
        xs = list(df['Fz (N) Vertical'])
        ys = list(df['Fx (N) Longitudinal'])

        return CubicSpline(xs, ys)

    # Creates a list of torque rpm pairs and their corresponding efficiency values to be interpolated from
    def createMotorData(self):
        df = pd.read_excel('AMK_Motor_Efficiency.xlsx', sheet_name="Sheet2")
        torque = list(df['Torque'])
        rpm = list(df['RPM'])
        eff = list(df['Efficiency'])
        points = []
        for i in range(len(torque)):
            points.append([torque[i], rpm[i]])
        self.motorPoints = points
        return eff

    # Returns motor efficiency from torque and rpm
    def interpolateEff(self, torque, rpm):
        return griddata(self.motorPoints, self.motorEff, (torque, rpm), method='cubic')

    # returns motor eff from aprx equation (less accurate but much faster)
    def effEq(self, torque, rpm):
        coeff = [0.5246, np.float32(8.303e-5), np.float32(0.08641), np.float32(-2.265e-8), np.float32(1.909e-5),
                 np.float32(-0.02905), np.float32(2.488e-12), np.float32(-1.748e-9), np.float32(-5.967e-7), np.float32(0.002848),
                 np.float32(-1.180e-16), np.float32(5.669e-14), np.float32(5.513e-11), np.float32(1.429e-9), np.float32(-1.254e-4),
                 np.float32(2.021e-21), np.float32(-4.058e-19), np.float32(-1.552e-15), np.float32(8.165e-14), np.float32(-9.192e-11), np.float32(2.076e-6)]
        count = 0
        index = 0
        eff = 0
        while count < 6:
            i = count
            while i > -1:
                eff += coeff[index] * (rpm ** i) * (torque ** (count - i))
                index += 1
                i -= 1
            count += 1

        if eff <= 0:
            eff = 0.01
        return eff

    # Converts torque and rpm to power used
    def motorPower(self, torque, rpm):
        eff = self.effEq(torque, rpm)
        return (2 * np.pi * torque * rpm / 60) / eff

    # Converts a power value into a corresponding torque request
    def powerToTorque(self, power, rpm):
        if rpm == 0:
            start = self.maxTorque
        else:
            start = (power * 60) / (4 * np.pi * rpm)
            if start > self.maxTorque:
                start = self.maxTorque

        t = np.arange(start, 0, -0.1)
        for i in t:
            if 2 * self.motorPower(i, rpm) < power:
                return i
        return 0

    # Takes 2 Cgs and their masses and outputs their combined CG
    def combineCG(self, cg1, mass1, cg2, mass2):
        cg = [0, 0, 0]

        for i in range(0, 3):
            cg[i] = (cg1[i] * mass1 + cg2[i] * mass2) / (mass1 + mass2)

        return cg

    # Finds the intersection of 2 lines from 3 sets of points per line
    # The first 2 points determine the slope of the line and the last point
    # gives a point of intersection for the line
    # Used in finding IC and PC
    def intersectFromPoints(self, ap1, ap2, ap3, bp1, bp2, bp3):
        x, z = symbols('x, z')
        a = self.createLine(ap1, ap2, ap3)
        b = self.createLine(bp1, bp2, bp3)
        output = solve([a, b], dict=True)
        return [output[0][x], 0, output[0][z]]

    # Creates where the first 2 points determine the slope of the line and
    # and passes through the third point
    # returns as a sympy equation
    def createLine(self, point1, point2, point3):
        x, z = symbols('x, z')
        slope = (point1[2] - point2[2]) / (point1[0] - point2[0])
        return sp.Eq(point3[2] + (slope * (x - point3[0])), z)
