import collections
import enum
import itertools
import time

import numpy as np
import PyMieScatt as ps
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

from pathlib import Path
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp


class Material:
    def __init__(self, k, rho, c, name="custom"):
        self.k = k
        self.rho = rho
        self.c = c
        self.rhoc = self.rho * self.c
        self.a = self.k / self.rhoc
        self.name = name

    def known():
        return {
            "Au": Material(k=317, rho=19300, c=129, name="Au"),
            "SiO2": Material(k=1.4, rho=2196, c=739, name="SiO2"),
            "H2O": Material(k=0.60652, rho=997, c=4181, name="H2O"),
            "air": Material(k=0.02638, rho=1.161, c=1007, name="air"),
        }

    def make(material):
        # Just return it if it is a Material
        if isinstance(material, Material):
            return material
        # Else, look it up in known_materials
        return Material.known()[material]


class Layer:
    def __init__(self, material, width, n_points):
        self.material = Material.make(material)
        self.width = width * 1e-9
        self.n_points = n_points


class Interface:
    def __init__(self, a, b, parent):
        self.kappa = a.material.k / b.material.k
        if a.material.name == b.material.name:
            self.resistance = 0
        else:
            self.resistance = parent.resistance_map[(a.material.name, b.material.name)]
        self.kapitza_length = (
            self.resistance * b.material.k / parent.unit_width
            + (a.dr / self.kappa + b.dr) / 2
        )


class Particle:
    def __init__(self, layers, resistance_map={}):
        self.layers = layers
        self.resistance_map = {
            ("Au", "SiO2"): 1 / 134e6,
            ("Au", "H2O"): 1 / 105e6,
            ("SiO2", "H2O"): 1 / 330e6,
            ("Au", "air"): 1 / 4e6,
            ("SiO2", "air"): 1 / 4e6,
        }
        self.resistance_map.update(resistance_map)
        self.initialize()

    def initialize(self):
        # keep track of location and number of first slice of each layer
        rr_0 = 0
        n_0 = 0
        # Units for normalizing
        self.unit_a = self.layers[-1].material.a
        self.unit_width = self.layers[0].width
        for layer in self.layers:
            layer.alpha = layer.material.a / self.unit_a
            layer.length = layer.width / self.unit_width
            layer.rr_0 = rr_0
            layer.dr = layer.length / layer.n_points
            layer.n_range = slice(n_0, n_0 + layer.n_points)
            # rr: locations of interfaces between slices
            layer.rr = np.linspace(rr_0, rr_0 + layer.length, layer.n_points + 1)
            layer.rr2 = layer.rr**2
            # r: locations of centers of slices, such that S(r)*dr = dV
            layer.r2 = (layer.rr[1:] ** 3 - layer.rr[:-1] ** 3) / (3 * layer.dr)
            layer.r = np.sqrt(layer.r2)
            layer.req = (layer.rr[1:] + layer.rr[:-1]) / 2
            layer.source = lambda t: 0
            layer.gamma = layer.alpha / layer.r2 / layer.dr
            rr_0 += layer.length
            n_0 += layer.n_points
        self.n_points = n_0
        self.r = np.concatenate([layer.r for layer in self.layers]) * self.unit_width
        self.req = (
            np.concatenate([layer.req for layer in self.layers]) * self.unit_width
        )
        self.interfaces = [
            Interface(a, b, self) for a, b in zip(self.layers, self.layers[1:])
        ]


Illumination = enum.Enum("Illumination", "NO FS NS CW")


class Parameters:
    def __init__(
        self,
        lambda0=532,
        isolated=False,
        time_end=1e-9,
        illumination=Illumination.FS,
        ode_options={},
        tau_ep0=1.7e-12,
        tau0=1e-9,
        irradiance=1e-3 / 1e-12,
        fluence=20,
        n_solvent=1.3330,
        cross_section=None,
    ):
        self.lambda0 = lambda0
        self.isolated = isolated
        self.illumination = illumination
        self.time_end = time_end
        self.ode_options = {
            "method": "BDF",
            "rtol": 1e-4,
            "atol": 1e-8,
            "dense_output": True,
        }
        self.ode_options.update(ode_options)
        self.tau_ep0 = tau_ep0
        self.tau0 = tau0
        self.irradiance = irradiance
        self.fluence = fluence
        self.n_solvent = n_solvent
        self.cross_section = cross_section


class Simulation:
    def __init__(self, particle, parameters):
        self.particle = particle
        self.parameters = parameters
        self.simulation_initialize()

    def simulation_initialize(self):
        core = self.particle.layers[0]
        self.time_unit = core.width**2 / self.particle.unit_a
        self.time_end = self.parameters.time_end / self.time_unit
        self.tau = self.parameters.tau0 / self.time_unit
        self.tau_ep = self.parameters.tau_ep0 / self.time_unit
        self.c_core = 4 / 3 * np.pi * core.width**3 * core.material.rhoc
        n_core = refractive_index(core.material.name)(self.parameters.lambda0)
        if self.parameters.cross_section is None:
            self.s_abs = (
                1e-18
                * ps.MieQ(
                    n_core,
                    self.parameters.lambda0,
                    2e9 * core.width,
                    self.parameters.n_solvent,
                    asCrossSection=True,
                )[2]
            )
        else:
            self.s_abs = self.parameters.cross_section
        # Set source in core according to Illumination type
        if self.parameters.illumination is Illumination.FS:
            core.source = (
                lambda t: self.time_unit
                * self.parameters.fluence
                * self.s_abs
                / self.c_core
                / self.parameters.tau_ep0
                * np.exp(-t / self.tau_ep)
            )
        elif self.parameters.illumination is Illumination.NS:
            core.source = (
                lambda t: self.time_unit
                * self.parameters.fluence
                * self.s_abs
                / self.c_core
                / self.parameters.tau0
                * np.exp(-np.pi * ((t - 2 * self.tau) ** 2 / self.tau**2))
            )
        elif self.parameters.illumination is Illumination.CW:
            core.source = (
                lambda t: self.time_unit
                * self.parameters.irradiance
                * self.s_abs
                / self.c_core
            )
        elif self.parameters.illumination is not Illumination.NO:
            raise NotImplementedError(
                f"Unknown illumination {self.parameters.illumination}"
            )

    def initialize(self):
        self.particle.initialize()
        self.simulation_initialize()

    def heat_eq(self, t, T):
        dTdt = np.zeros_like(T)
        Tslices = [np.zeros(layer.n_points) for layer in self.particle.layers]
        dTdr = [np.zeros(layer.n_points + 1) for layer in self.particle.layers]
        for i, layer in enumerate(self.particle.layers):
            Tslices[i] = T[layer.n_range]
            dTdr[i][1:-1] = np.diff(Tslices[i]) / layer.dr
        for i, interface in enumerate(self.particle.interfaces):
            dTdr[i + 1][0] = (
                Tslices[i + 1][0] - Tslices[i][-1]
            ) / interface.kapitza_length
            dTdr[i][-1] = dTdr[i + 1][0] / interface.kappa
        if not self.parameters.isolated:
            dTdr[-1][-1] = -T[-1] / self.particle.layers[-1].dr
        for layer, Tr in zip(self.particle.layers, dTdr):
            dTdt[layer.n_range] = layer.gamma * np.diff(layer.rr2 * Tr) + layer.source(
                t
            )
        return dTdt

    def run(self):
        T = np.zeros(self.particle.n_points)
        start = time.perf_counter()
        solution = solve_ivp(
            self.heat_eq, (0, self.time_end), T, **self.parameters.ode_options
        )
        r = self.particle.req * 1e9
        sol = lambda t: xr.DataArray(
            solution.sol(np.asarray(t) / self.time_unit / 1e12),
            [("space", r), ("time", t)],
            name="temperature (°C)",
        )
        solve_time = time.perf_counter() - start
        attrs = {
            "simulation": self,
            "sol": sol,
            "solution": solution,
            "solve_time": solve_time,
            "particle": self.particle,
            "parameters": self.parameters,
        }
        result = xr.DataArray(
            solution.y,
            [("space", r), ("time", solution.t * self.time_unit * 1e12)],
            name="temperature (°C)",
            attrs=attrs,
        )
        return result


def refractive_index(material_name):
    conversion = 1239.841984
    n = np.loadtxt(Path(__file__).with_name(f"n{material_name}.txt"))
    inter = interp1d(n[:, 0], np.log(n[:, 1:]), "cubic", 0)
    return lambda x: complex(*np.exp(inter(conversion / x)))
