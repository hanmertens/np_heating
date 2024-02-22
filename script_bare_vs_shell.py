#!/usr/bin/env python3

from np_heating import *

# Nanoparticle core radius in nanometer (e.g. 15 or 30)
radius = 15
# Medium (choose from H2O, air, EtOH, toluene, hexane or CHCl_3)
medium = "H2O"
# Simulate long timescale (True: 100 ns; False: 800 ps)
long_timescale = False

gold = "#ffb81c"
gray = "#808080"


def get_results(r, mater, resistance_map={}, **kwargs):
    n = r // 15 * 20
    if mater == "air":
        particle_b = Particle(
            [
                Layer("Au", r, n),
                Layer("air", 20, 200),
                Layer("air", 80, 200),
                Layer("air", 1600, 1000),
            ],
            resistance_map=resistance_map,
        )
        particle_s = Particle(
            [
                Layer("Au", r, n),
                Layer("SiO2", 2, 40),
                Layer("air", 20, 200),
                Layer("air", 80, 200),
                Layer("air", 1600, 1000),
            ],
            resistance_map=resistance_map,
        )
    else:
        particle_b = Particle(
            [
                Layer("Au", 15, 20),
                Layer(mater, 20, 200),
                Layer(mater, 80, 200),
                Layer(mater, 160, 100),
            ],
            resistance_map=resistance_map,
        )
        particle_s = Particle(
            [
                Layer("Au", 15, 20),
                Layer("SiO2", 2, 40),
                Layer(mater, 20, 200),
                Layer(mater, 80, 200),
                Layer(mater, 160, 100),
            ],
            resistance_map=resistance_map,
        )
    param_kwargs = dict(time_end=0.8e-9, fluence=16, lambda0=400)
    param_kwargs.update(kwargs)
    parameters = Parameters(**param_kwargs)
    results = {}
    simulation_b = Simulation(particle_b, parameters)
    simulation_s = Simulation(particle_s, parameters)
    results["b"] = simulation_b.run()
    results["s"] = simulation_s.run()
    nrange = particle_s.layers[1].n_range
    results["nrange"] = nrange
    results["ninter"] = [nrange.start - 1, nrange.start, nrange.stop - 1, nrange.stop]
    return results


def plot_temporal(result, shell, medium, ninter):
    medium = medium.title()
    if shell:
        result.isel(space=ninter).plot.line(x="time")
        layers = ["Core", "Shell (inner)", "Shell (outer)", medium]
    else:
        result.isel(space=ninter[0]).plot.line(x="time", color="tab:blue")
        result.isel(space=ninter[1]).plot.line(x="time", color="tab:red")
        plt.title("")
        layers = ["Core", medium]
    plt.legend(layers, title="Interface")
    plt.ylabel("$\Delta T$ ($\degree$C)")


def plot_radial(result, taxis, shell, cmap, markerdict):
    plt.gca().set_prop_cycle(color=cmap)
    data = result.sol(taxis)
    particle = result.particle
    for layer in particle.layers:
        n_range = layer.n_range
        data.isel(space=n_range).plot.line(x="space")
    plt.legend([f"{t}" for t in taxis], title="$t$ (ps)", ncols=2)
    markery = plt.ylim()[1] * 1.025
    corew = result.particle.layers[0].width * 1e9
    mediumw = corew
    plt.axvspan(0, corew, color=gold, alpha=0.3, linewidth=0)
    if shell:
        shellw = result.particle.layers[1].width * 1e9 + corew
        mediumw = shellw
        plt.axvspan(corew, shellw, color=gray, alpha=0.3, linewidth=0)
        plt.plot(corew, markery, fillstyle="right", color="tab:orange", **markerdict)
        plt.plot(shellw, markery, fillstyle="left", color="tab:green", **markerdict)
    plt.plot(corew, markery, fillstyle="left", color="tab:blue", **markerdict)
    plt.plot(mediumw, markery, fillstyle="right", color="tab:red", **markerdict)


def axis_off(which):
    for i in which:
        dict(x=plt.xlabel, y=plt.ylabel)[i]("")
        plt.tick_params(i, labelleft=False, labelbottom=False)


def rlabel(label):
    plt.gca().yaxis.set_label_position("right")
    plt.ylabel(label)


def plot_all(results, tlim, rlim, medium, taxis, cmap):
    markerdict = dict(marker="v", markeredgewidth=0, markersize=10, clip_on=False)
    ninter = results["ninter"]

    tplot = plt.subplot(221)
    plot_temporal(results["b"], False, medium, ninter)
    plt.xlim(tlim)
    plt.ylim(bottom=0)
    axis_off("x")
    plt.title("(a)", loc="left", weight="bold")

    plt.subplot(223, sharex=tplot, sharey=tplot)
    plot_temporal(results["s"], True, medium, ninter)
    plt.title("(c)", loc="left", weight="bold")

    rplot = plt.subplot(222, sharey=tplot)
    plot_radial(results["b"], taxis, False, cmap, markerdict)
    plt.xlim(rlim)
    axis_off("xy")
    rlabel(f"Bare in {medium}")
    plt.title("(b)", loc="left", weight="bold")

    plt.subplot(224, sharex=rplot, sharey=tplot)
    plot_radial(results["s"], taxis, True, cmap, markerdict)
    axis_off("y")
    rlabel(f"Core$-$shell in {medium}")
    plt.title("(d)", loc="left", weight="bold")


taxis = [10, 20, 40, 100, 200, 400]
cmap = ["tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]

etoh = Material(k=0.167, rho=789.3, c=2438, name="EtOH")
toluene = Material(k=0.1310, rho=862.3, c=1707, name="toluene")
hexane = Material(k=0.1167, rho=660.6, c=2270, name="hexane")
chcl3 = Material(k=0.117, rho=1478.8, c=957, name="CHCl_3")
resistance_map = {
    ("Au", "EtOH"): 1 / 32e6,
    ("SiO2", "EtOH"): 1 / 87e6,
    ("Au", "toluene"): 1 / 30e6,
    ("SiO2", "toluene"): 1 / 87e6,
    ("Au", "hexane"): 1 / 25e6,
    ("SiO2", "hexane"): 1 / 87e6,
    ("Au", "CHCl_3"): 1 / 30e6,
    ("SiO2", "CHCl_3"): 1 / 87e6,
}

print(f"NP core radius: {radius} nm")
print(f"Medium: {medium}")
print(f"Long timescale: {long_timescale}")

medium_dict = dict(
    EtOH=Material(k=0.167, rho=789.3, c=2438, name="EtOH"),
    toluene=Material(k=0.1310, rho=862.3, c=1707, name="toluene"),
    hexane=Material(k=0.1167, rho=660.6, c=2270, name="hexane"),
    CHCl_3=Material(k=0.117, rho=1478.8, c=957, name="CHCl_3"),
)
medium = medium_dict.get(medium, medium)

if long_timescale:
    kwargs = dict(time_end=100e-9)
    tlim = [0, 10_000]
else:
    kwargs = {}
    tlim = [0, 400]

results = get_results(radius, medium, resistance_map, **kwargs)

s_abs = 1e18 * results["b"].simulation.s_abs
print(f"NP core absorbance cross section: {s_abs} nm^2")

plt.figure(figsize=(8, 6), dpi=300)
plot_all(results, tlim, [0, 50], medium, taxis, cmap)
plt.savefig("script_bare_vs_shell.png")
