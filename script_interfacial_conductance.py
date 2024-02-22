#!/usr/bin/env python3

from np_heating import *

# Nanoparticle core radius in nanometer (e.g. 15 or 30)
radius = 15
# Shell thickness in nanometer (e.g. 0 or 4)
thickness = 0

print(f"NP core radius: {radius} nm")
if thickness == 0:
    print(f"NP shell thickness: {thickness} nm (no shell)")
else:
    print(f"NP shell thickness: {thickness} nm")

layers = [Layer("Au", 15, 20)]
if thickness > 0:
    layers.append(Layer("SiO2", thickness, 10 * thickness))
layers += [
    Layer("H2O", 20, 200),
    Layer("H2O", 80, 200),
    Layer("H2O", 1600, 1000),
]
if thickness == 0:
    gees = [45, 105, 250, 450]
    rmaps = [{("Au", "H2O"): 1 / (g * 1e6)} for g in gees]
else:
    gees = [(39, 220), (39, 2080), (222, 220), (222, 2080)]
    rmaps = [
        {("Au", "SiO2"): 1 / (g[0] * 1e6), ("SiO2", "H2O"): 1 / (g[1] * 1e6)}
        for g in gees
    ]
particles = [Particle(layers, rmap) for rmap in rmaps]
parameters = Parameters(time_end=0.8e-9, fluence=10, lambda0=400)
simulations = [Simulation(particle, parameters) for particle in particles]
results = [simulation.run() for simulation in simulations]

nrange = particles[0].layers[1].n_range
ninter = [nrange.start - 1, nrange.start, nrange.stop - 1, nrange.stop]

plt.figure(figsize=(8, 6), dpi=300)
for i, result in enumerate(results):
    plt.subplot(221 + i)
    if thickness == 0:
        result.isel(space=ninter[0]).plot.line(x="time", color="tab:blue")
        result.isel(space=ninter[1]).plot.line(x="time", color="tab:red")
        plt.legend(["Core", "Water"], loc="upper right")
        plt.title(f"$G$ = {gees[i]} MW m$^{{-2}}$ K$^{{-1}}$")
    else:
        result.isel(space=ninter).plot.line(x="time")
        plt.legend(["Core", "Shell (inner)", "Shell (outer)", "Water"])
        plt.title(f"$G$ = {gees[i]} MW m$^{{-2}}$ K$^{{-1}}$")
    plt.xlim(0, 400)
    plt.ylim(0, 250)
    plt.xlabel("$t$ (ps)")
    plt.ylabel("$\Delta T$ ($\degree$C)")
    plt.title(f"({chr(ord('a')+i)})", loc="left", weight="bold")
    if i % 2:
        plt.ylabel("")
        plt.tick_params("y", labelleft=False)
    if i < 2:
        plt.xlabel("")
        plt.tick_params("x", labelbottom=False)
plt.savefig("script_interfacial_conductance.png")
