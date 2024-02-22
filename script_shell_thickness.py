#!/usr/bin/env python3

from np_heating import *

# Nanoparticle core radius in nanometer (e.g. 15 or 30)
radius = 15
# Shell thicknesses in nanometer (e.g. a np.arange from 0 to 16)
thickness = np.arange(16)

print(f"NP core radius: {radius} nm")
print(f"NP shell thicknesses: {thickness} nm")


def coreshell(n):
    """Create Au@SiO2 NP with a shell thickness of n nm."""
    base = [Layer("Au", radius, radius // 15 * 20)]
    if n > 0:
        base.append(Layer("SiO2", n, int(10 * n)))
    base += [Layer("H2O", 20, 200), Layer("H2O", 80, 200), Layer("H2O", 160, 100)]
    particle = Particle(base)
    return particle


def np_invarea(d):
    return 1 / (4 * np.pi * (15 + d) ** 2)


parameters = Parameters(time_end=0.8e-9, fluence=10, lambda0=400)
interfaces = 20 + 10 * thickness
full_result = []
for n in thickness:
    particle = coreshell(n)
    simulation = Simulation(particle, parameters)
    thisresult = simulation.run()
    full_result.append(thisresult)

t = np.arange(800, step=0.1)
data = np.zeros((len(thickness), len(t)))
for i, (ri, res) in enumerate(zip(interfaces, full_result)):
    data[i, :] = res.sol(t).isel(space=ri)
tdat = xr.DataArray(data, [("thickness", thickness), ("time", t)])

fraction = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
time = xr.DataArray(
    np.zeros([len(thickness), len(fraction)]),
    [("thickness", thickness), ("fraction", fraction)],
)
temp = xr.DataArray(
    np.zeros([len(thickness), len(fraction)]),
    [("thickness", thickness), ("fraction", fraction)],
)
for r in thickness:
    rdata = tdat.sel(thickness=r)
    rtemp = rdata.max()
    for f in fraction:
        fdata = rdata.where(rdata >= rtemp * f, drop=True).isel(time=0)
        time.loc[r, f] = fdata.time
        temp.loc[r, f] = fdata.data

stime = time.where(time.thickness > 0, drop=True)
timefit = stime.polyfit(dim="thickness", deg=1)
stemp = temp.where(temp.thickness, drop=True).where(temp.fraction > 0, drop=True)
orithickness = stemp.thickness.data
stemp["thickness"] = np_invarea(orithickness)
stemp = stemp.assign_coords(orithickness=("thickness", orithickness))
tempfit = stemp.polyfit(dim="thickness", deg=1)
detthickness = np.linspace(0.5, np.max(thickness) + 0.5, 100)
detthickness = xr.DataArray(
    detthickness, [("thickness", detthickness)], name="thickness"
)

totfraction = np.asarray(100 * np.asarray(fraction), dtype=int)
cmap = [mpl.cm.viridis(i / (len(thickness) - 1)) for i in range(len(thickness))]
cmap2 = [mpl.cm.plasma(i / (len(totfraction) - 1)) for i in range(len(totfraction))]

plt.figure(figsize=(8, 6), dpi=300)
fig1 = plt.subplot(211)
plt.gca().set_prop_cycle(color=cmap)
tdat.plot.line(x="time")
plt.xlim(0, 400)
plt.ylim(bottom=0)
leg1 = plt.legend(
    thickness,
    title="Shell thickness (nm)",
    ncol=4,
    bbox_to_anchor=(0.5, 1),
    loc="lower right",
)
plt.gca().set_prop_cycle(color=cmap2)
for f in range(len(totfraction)):
    marker = "o" if f == len(totfraction) - 1 else "."
    plt.plot(time[:, f], temp[:, f], marker=marker, lw=0)
plt.xlabel("$t$ (ps)")
plt.ylabel("$\Delta T$ ($\degree$C)")
plt.title("(a)", loc="left", weight="bold")

fig2 = plt.subplot(224, sharey=fig1)
plt.gca().set_prop_cycle(color=cmap2[-1:])
plt.plot(temp.sel(fraction=1), "o")
plt.plot(
    detthickness,
    xr.polyval(np_invarea(detthickness), tempfit.polyfit_coefficients).sel(fraction=1),
    lw=0.5,
)
plt.xlabel("Shell thickness (nm)")
plt.ylabel("$\Delta T_\mathrm{max}$ ($\degree$C)")
plt.title("(c)", loc="left", weight="bold")

fig3 = plt.subplot(223, sharex=fig2)
plt.gca().set_prop_cycle(color=cmap2)
dots = plt.plot(time[:, :-1], ".")
circles = plt.plot(time[:, -1], "o")
plt.plot(detthickness, xr.polyval(detthickness, timefit.polyfit_coefficients), lw=0.5)
plt.xlabel("Shell thickness (nm)")
plt.ylabel("$t$ (ps)")
plt.xlim([-0.5, 15.5])
plt.xticks([0, 5, 10, 15])
plt.ylim(bottom=0)
fig1.legend(
    dots + circles,
    totfraction,
    title="$\Delta T$ (% of $\Delta T_\mathrm{max}$)",
    bbox_to_anchor=(0.5, 1),
    loc="lower left",
    ncol=3,
)
fig1.add_artist(leg1)
plt.title("(b)", loc="left", weight="bold")
plt.savefig("script_shell_thickness.png")
