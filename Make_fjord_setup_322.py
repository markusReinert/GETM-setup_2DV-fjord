"""Create configuration 322 of the 2D model for the 79NG fjord.

Fjord configuration 322:
 – 3: hyperbolic tangent lower ice edge
 – 2: parabolic sill
 – 2: deep trough
(Earlier model configurations had other shapes for the topography,
which were labeled with lower numbers.)

The bottom topography consists of three or four parts, which are
from left right, i.e., from the grounding line to the open ocean:
(1) A deep trough that begins at the grounding line at a given depth
    with zero slope, reaches a given maximum depth at a given point,
    and ends where its derivative reaches a given maximum slope.
(2) A linear part with the given (maximum) slope.
(3) An optional parabolic sill that reaches its highest point
    at a given minimum depth at a given position.
(4) An exponential continental shelf seaward of the parabolic sill
    that approaches a given depth far offshore.
These parts are connected in a way that the derivative is continuous.

The ice topography has a hyperbolic tangent shape with a given maximum
slope at the grounding line.
After the calving front, the ice goes linearly to zero.

The 3D boundary conditions (i.e. stratification) can be time-dependent.

The 2D boundary condition is a constant zero elevation.

All x,y,z-values in this script are given in meter.

by Markus Reinert and Sarah Vierling at IOW, November 2020 to March 2023
"""

import os
from datetime import datetime
import xml.etree.ElementTree as ET

import numpy as np
import xarray as xr
from scipy.optimize import fsolve
from lib.eos import rho, T_freezing


OUTPUT_PATH = "."


XML_file = "fjord_322.xml"
print(f"Loading GETM settings from {XML_file}.")
XML_tree = ET.parse(XML_file)
XML_root = XML_tree.getroot()
assert XML_root.tag == "scenario", 'root-element of XML-file is not "scenario".'

g = float(XML_root.find("getm/parameters/g").text)
rho_ice = float(XML_root.find("getm/ice/rho_ice").text)
dt = float(XML_root.find("getm/time/timestep").text)


# Define horizontal grid
# Set resolution in along-fjord direction
dx = 500.0
# Set width in across-fjord direction (cf. Mayer et al. 2000)
dy = 20e3
# Define x-axis at center points with first point a land point and water for x >= 0
xax = np.arange(-dx/2, 150e3, dx)
# Define y-axis at center points with one water point from 0 to dy, surrounded by two land points
yax = np.array([-dy/2, dy/2, 3*dy/2])

print(f"Grid size: {xax.size} × {yax.size}")
print(f"Grid resolution: {dx} m × {dy/1e3} km")

# Create the output arrays for the fjord configuration
elev = np.zeros((yax.size, xax.size))
glIceD = np.zeros((yax.size, xax.size))
bathymetry = np.zeros((yax.size, xax.size))

# Mark boundary points (i.e. land points)
bathymetry[0, :] = np.nan
bathymetry[1, 0] = np.nan
bathymetry[2, :] = np.nan


# Define bottom topography

# Grounding line
XML_gline = XML_root.find("fjord/bathymetry/grounding_line")
# x_gline = 0 (used only implicitly)
z_gline = -float(XML_gline.get("depth"))
# Deepest point of the fjord in the trough
XML_trough = XML_root.find("fjord/bathymetry/trough")
x_min = float(XML_trough.get("x"))
z_min = -float(XML_trough.get("depth"))
# Highest point of the fjord on the sill
XML_sill = XML_root.find("fjord/bathymetry/sill")
x_max = float(XML_sill.get("x"))
z_max = -float(XML_sill.get("depth"))
# Continental shell seaward of the sill
XML_shelf = XML_root.find("fjord/bathymetry/shelf")
# x_shelf -> +infinity (theoretically)
z_shelf = -float(XML_shelf.get("depth"))
# Maximum slope of the bathymetry between the deepest and highest point
dzdx_max = float(XML_root.find("fjord/bathymetry/max_slope").text)

# We construct a deep trough bathymetry of cubic form: z = ax³+bx²+cx+d.
# We require that the bottom slope at the grounding line (x=0) is zero,
# so c=0, thus: z = ax³+bx²+d.
# We prescribe the depth of the grounding line, so d=z_gline, thus: z = ax³+bx²+z_gline
# From the condition that the bottom slope at the minimum is zero, we get:
b = (3/(x_min**2))*(z_min - z_gline)
a = (-2/3)*(b/x_min)
trough = lambda x: a * x**3 + b * x**2 + z_gline

# From the condition that the deep trough ends where the slope equals
# dzdx_max, we get its seaward end point:
x0 = ((-2*b) + np.sqrt((2*b)**2 + 12*a*dzdx_max)) / (6*a)
z0 = trough(x0)

# Define linear connection between trough and sill/shelf
slope = lambda x: z0 + dzdx_max * (x - x0)

bathymetry[1, xax > 0]  = trough(xax[xax > 0])
bathymetry[1, xax > x0] = slope(xax[xax > x0])

# When the sill has the same depth as the shelf, there is no sill
if z_max == z_shelf:
    # Neglect the sill
    z1 = -600  # z-limit of linear connection (i.e. where shelf begins)
    x1 = (z1 - (z0 - (dzdx_max*x0))) / dzdx_max  # x-limit of linear connection
    b2 = dzdx_max / (z_shelf - z1)
    a2 = (z_shelf - z1) / (np.exp(-b2*x1))
    shelf = lambda x: -a2*np.exp(-x*b2) + z_shelf

    bathymetry[1, xax > x1] = shelf(xax[xax > x1])

else:
    # Incorporate the sill
    # Calculate the landward end point of the sill
    x1 = 2*(z_max - z0) / dzdx_max + 2*x0 - x_max
    z1 = dzdx_max*x1 + (z0 - (dzdx_max*x0))
    # Calculate the curvature of the sill
    a_sill = dzdx_max / (x1 - x_max)
    sill = lambda x: a_sill / 2 * (x - x_max)**2 + z_max

    # Define seaward end point of the sill
    if z_max <= -400:
        # … such that the slope is one third of the maximum slope
        x2 = x_max - 1/3 * dzdx_max / a_sill
    else:
        # … such that the slope is half of the maximum slope
        x2 = x_max - 1/2 * dzdx_max / a_sill
    z2 = sill(x2)

    print("Parabolic sill")
    print(f"  starts at  ({x1:_.0f}, {z1:.1f}),")
    print(f"  is max. at ({x_max:_.0f}, {z_max:.1f}),")
    print(f"  ends at    ({x2:_.0f}, {z2:.1f}),")
    print(f"  and covers {(x2-x1)/dx:.1f} grid points.")

    # Define an exponential shelf seaward of the sill
    b_shelf = a_sill * (x2 - x_max) / (z2 - z_shelf)
    a_shelf = (z2 - z_shelf) / np.exp(b_shelf * x2)
    shelf = lambda x: a_shelf * np.exp(b_shelf * x) + z_shelf
    
    bathymetry[1, xax > x1] = sill(xax[xax > x1])
    bathymetry[1, xax > x2] = shelf(xax[xax > x2])

# Make shelf constant near the boundary
bathymetry[1, -10:] = bathymetry[1, -10] 


print("Deepest point of the bathymetry: ", end="")
print(f"({xax[np.nanargmin(bathymetry[1])]:_.0f}, {-np.nanmin(bathymetry):.1f})")

print("Largest possible timestep according to CFL criterion in 2D (ignoring ice):")
c = np.sqrt(g * -np.nanmin(bathymetry))
dt_max = dx / c 
print(f"  dt_max = dx / c = ({dx} m) / ({c:.2f} m/s) = {dt_max:.2f} s")
assert dt < dt_max, "CFL criterion violated, GETM timestep too large!"

print("Slope of the bathymetry")
print(f"  should be at most {dzdx_max} in absolute value, ", end="")
print(f"is between {np.min(np.diff(bathymetry[1, 1:]) / dx):.4f} ", end="")
print(f"and {np.max(np.diff(bathymetry[1, 1:]) / dx):.4f}.")


# Define ice topography: hyperbolic tangent shape z = a*tanh(b*(x-c))+d,
# with a,b,c,d defined by the depth at the grounding line, the depth far
# offshore, the maximum slope, and the x-position of the maximum slope.

# Ice depth at the grounding line
z_ice_gline = z_gline
XML_calving = XML_root.find("fjord/ice_tongue/calving_front")
# Position of the calving front; comparing the map of Schaffer et al.
# (2020) with the distance along the section by Mayer et al. (2000),
# the length of the floating ice tongue should be between 70 and 75 km;
# BedMachine data suggest rather 70 km or even less.
x_ice_calving = float(XML_calving.get("x"))
# Ice depth at the calving front; Schaffer et al. (2020) show it as
# about 100 m, Mayer et al. (2000) show it less deep; looking at
# BedMachine data, there are even several places with only 50 m depth,
# so it should be between 50 and 100 m
z_ice_calving = -float(XML_calving.get("depth"))
# Slope of the calving front; in reality, the calving front is almost vertical,
# but this is not possible with terrain-following coordinates,
# so we replace the vertical front by a linear front
dzdx_calving = float(XML_root.find("fjord/ice_tongue/calving_front_slope").text)
# Maximum slope of the ice shape
dzdx_max_ice = float(XML_root.find("fjord/ice_tongue/max_slope").text)
# Position of the maximum slope (e.g.: 10e3 for 10 km from grounding line)
c_ice = 0

if c_ice == 0:
    # Calculate all parameters analytically
    a_ice = z_ice_calving - z_ice_gline
else:
    # Calculate the parameter a numerically, the other parameters analytically
    a_ice = fsolve(
        lambda a: a * (1 + np.tanh(dzdx_max_ice / a * c)) - z_ice_calving + z_ice_gline,
        z_ice_calving - z_ice_gline,
    )
b_ice = dzdx_max_ice / a_ice
d_ice = z_ice_gline + a_ice * np.tanh(b_ice * c_ice)

# Define the ice shape up to the calving front with the desired ice shape
elev[1][xax <= x_ice_calving] = (
    a_ice * np.tanh(b_ice * (xax[xax <= x_ice_calving] - c_ice)) + d_ice
)
# Add a linear slope from the calving front to the ocean surface
elev[1][xax > x_ice_calving] = (
    z_ice_calving + dzdx_calving * (xax[xax > x_ice_calving] - x_ice_calving)
)
# Replace positive values with zero
elev[1][elev[1] > 0] = 0

print("Slope of the ice topography")
print(f"  should be at most {dzdx_max_ice}, is at most {np.max(np.diff(elev[1]) / dx):.4f}.")


# The initial model stratification is defined by layers at a given
# depth (-z) with given values of temperature (T) and salinity (S).
# Between these z-levels, T and S are linearly interpolated.
# The vertical boundary condition is equal to the initial stratification.
# An option for time-varying boundary conditions is implemented below.
stratification = []
XML_stratification = XML_root.find("fjord/stratification")
for level in XML_stratification:
    assert level.tag == "level",\
        f"all elements in stratification must be named 'level', not {child.tag!r}"
    z_level = float(level.get("z"))
    S_level = float(level.get("S"))
    T_level = float(level.get("T"))
    stratification.append((z_level, S_level, T_level))
stratification.sort(key=lambda level: level[0])
z_discrete = [level[0] for level in stratification]
S_discrete = [level[1] for level in stratification]
T_discrete = [level[2] for level in stratification]


assert min(z_discrete) <= np.nanmin(bathymetry), "bathymetry deeper than stratification"


# Specify how the open ocean stratification changes with time, or set to
# empty lists for constant-in-time boundary conditions.  If the lists
# are non-empty, they all must have the same size.  The profiles of T
# and S at the open boundary (i.e., the open ocean) stay constant until
# the first datetime value in time_BC, and they stay constant after the
# last.  In between, they change linearly.  Each value of T_in_time and
# S_in_time must have the same size as z_discrete.
time_BC = []    # [datetime(2000, 2, 1), datetime(2000, 5, 1)]
T_in_time = []  # [T_discrete, T_discrete]
S_in_time = []  # [S_discrete, [35, 35, 34]]


# Compute ice thickness
if (
        all(T_value == T_discrete[0] for T_value in T_discrete) and
        all(S_value == S_discrete[0] for S_value in S_discrete) and
        len(time_BC) == 0
):
    rho_value = rho(S_discrete[0], T_discrete[0], p=0)
    print("Unstratified model with density {:.2f} kg/m³".format(rho_value))
    glIceD = rho_value * -elev / rho_ice
    # Check that temperature is above freezing point
    assert T_discrete[0] >= T_freezing(S_discrete[0], 0),\
        "water temperature below freezing point at the surface"
    # In the unstratified case, the calculation of the glacial ice
    # thickness is exact, which corresponds to dz -> 0
    dz = 0.0
else:
    # Set accuracy of calculations involving density
    dz = 0.001
    # Define a vertical axis for the integration
    z_profile = np.arange(-dz/2, np.nanmin(bathymetry), -dz)
    print(f"z-axis goes from {z_profile.min()} m to {z_profile.max()} m in steps of {dz} m.")

    # Calculate initial profiles of salinity, temperature, and density
    S_profile = np.interp(z_profile, z_discrete, S_discrete)
    T_profile = np.interp(z_profile, z_discrete, T_discrete)
    rho_profile = rho(S_profile, T_profile, p=0)

    # Check that initial temperature profile is nowhere less than freezing point
    assert np.all(T_profile >= T_freezing(S_profile, z_profile)),\
        "water temperature below local freezing point"
    # Check for vertical stability
    assert np.all(np.diff(rho_profile) >= 0), "initial stratification not stable"
    # Idea: these two checks could also be made for time-varying BCs

    # Compute thickness of glacial ice: using Archimedes' principle, we
    # calculate how much ice is needed to keep the ocean surface at the
    # initial surface elevation
    for j in range(yax.size):
        for i in range(xax.size):
            glIceD[j, i] = np.sum(rho_profile[z_profile >= elev[j, i]] * dz) / rho_ice


# Create and save dataset with fjord configuration
ds = xr.Dataset(
    {
        "grid_type": 1,
        "dx": ([], dx, {"long_name": "grid spacing in x-direction", "units": "m"}),
        "dy": ([], dy, {"long_name": "grid spacing in y-direction", "units": "m"}),
        "bathymetry": (
            ["yax", "xax"],
            -bathymetry,
            {"long_name": "bottom depth", "units": "m", "positive": "down"},
        ),
        "elev": (
            ["yax", "xax"],
            elev,
            {"long_name": "initial surface elevation", "units": "m", "positive": "up"},
        ),
        "glIceD": (
            ["yax", "xax"],
            glIceD,
            {"long_name": "glacial ice thickness", "units": "m"},
        ),
    },
    coords={
        "xax": (["xax"], xax, {"long_name": "x-axis", "units": "m"}),
        "yax": (["yax"], yax, {"long_name": "y-axis", "units": "m"}),
    },
    attrs={
        "title": "Configuration 322 of the 79NG fjord model for GETM",
        "author": "Markus Reinert & Sarah Vierling (IOW)",
        "rho_ice": float(rho_ice),
        "T-stratification": ", ".join([
            f"{T_val} degC at {z_val} m" for T_val, z_val in zip(T_discrete, z_discrete)
        ]),
        "S-stratification": ", ".join([
            f"{S_val} psu at {z_val} m" for S_val, z_val in zip(S_discrete, z_discrete)
        ]),
        "dz": f"{dz} m (accuracy used in the computation of glIceD)",
    },
)

filename = os.path.join(OUTPUT_PATH, "fjord_322.nc")
ds.to_netcdf(
    filename,
    encoding={
        "dx": {"_FillValue": None},
        "dy": {"_FillValue": None},
        "xax": {"_FillValue": None},
        "yax": {"_FillValue": None},
        "elev": {"_FillValue": None},
        "glIceD": {"_FillValue": None},
        "bathymetry": {"_FillValue": None},
    }
)
print(f"Saved fjord configuration as {filename!r}.")


# Create and save dataset with 3D boundary conditions
assert len(time_BC) == len(S_in_time) == len(T_in_time)
if len(time_BC) == 0:
    S_matrix = np.ones((2, 1, len(z_discrete))) * S_discrete[::-1]
    T_matrix = np.ones((2, 1, len(z_discrete))) * T_discrete[::-1]
    time_list = [datetime(2000, 1, 1), datetime(2100, 1, 1)]
else:
    S_matrix = [
        [S_discrete[::-1]],
        *([S[::-1]] for S in S_in_time),
        [S_in_time[-1][::-1]],
    ]
    T_matrix = [
        [T_discrete[::-1]],
        *([T[::-1]] for T in T_in_time),
        [T_in_time[-1][::-1]],
    ]
    time_list = [
        datetime(2000, 1, 1),
        *time_BC,
        datetime(2100, 1, 1),
    ]
ds = xr.Dataset(
    {
        "salt": (
            ["time", "nbdy", "zax"],
            S_matrix,
            {"long_name": "salinity", "units": "PSU"},
        ),
        "temp": (
            ["time", "nbdy", "zax"],
            T_matrix,
            {"long_name": "temperature", "units": "degC"},
        ),
    },
    coords={
        "nbdy": (["nbdy"], [0]),
        "time": (["time"], time_list),
        "zax": (
            ["zax"],
            z_discrete[::-1],
            {"long_name": "z-axis", "units": "m", "positive": "up"},
        ),
    },
    attrs={
        "title": "Boundary conditions (3D) for the 79NG fjord",
        "author": "Markus Reinert (IOW)",
        "description": "Boundary conditions are given by piece-wise linear z-profiles of S and T.",
    },
)
filename = os.path.join(OUTPUT_PATH, "bdy_3d.nc")
ds.to_netcdf(
    filename,
    unlimited_dims=["time"],
    encoding={
        "salt": {"_FillValue": None},
        "temp": {"_FillValue": None},
        "time": {"units": "seconds since 2000-01-01"},
    },
)
print(f"Saved 3D boundary conditions as {filename!r}.")


# Create and save dataset with 2D boundary conditions
ds = xr.Dataset(
    {
        "elev": (
            ["time", "nbdy"],
            np.zeros((2, 1)),
            {"long_name": "sea surface elevation", "units": "m"},
        ),
    },
    coords={
        "nbdy": (["nbdy"], [0]),
        "time": (["time"], [datetime(2000, 1, 1), datetime(2100, 1, 1)]),
    },
    attrs={
        "title": "Boundary conditions (2D) for the 79NG fjord",
        "author": "Markus Reinert (IOW)",
        "description": "Boundary conditions are given by constant zero elevation.",
    },
)
filename = os.path.join(OUTPUT_PATH, "bdy_2d.nc")
ds.to_netcdf(
    filename,
    unlimited_dims=["time"],
    encoding={
        "elev": {"_FillValue": None},
        "time": {"units": "seconds since 2000-01-01"},
    },
)
print(f"Saved 2D boundary conditions as {filename!r}.")


# Write file with salinity profile
filename = os.path.join(OUTPUT_PATH, "salt_profile.txt")
with open(filename, "w") as f:
    f.write(f"{len(z_discrete)}\n")
    for z_val, S_val in zip(reversed(z_discrete), reversed(S_discrete)):
        f.write(f"{z_val:6} {S_val:3}\n")
print(f"Saved initial salinity stratification as {filename!r}.")


# Write file with temperature profile
filename = os.path.join(OUTPUT_PATH, "temp_profile.txt")
with open(filename, "w") as f:
    f.write(f"{len(z_discrete)}\n")
    for z_val, T_val in zip(reversed(z_discrete), reversed(T_discrete)):
        f.write(f"{z_val:6} {T_val:5}\n")
print(f"Saved initial temperature stratification as {filename!r}.")
