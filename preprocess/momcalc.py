"""
Created on Thu Aug  1 14:39:10 2019

@author: zhangyurong
"""

import os
import numpy as np
import obspy.taup


def correct(dist, z_src, phase):
    """
    Inputs:
        - mod_name = name of velocity model (e.g. iasp91)
        - dist = epicentral distance (degrees)
        - z_src = source depth (km)
        - phase = seismic phase (currently supports body waves "P" and "S";
          needs to be upper case)
    Internal parameters:
        - v[p|s]_[s|r] = P-wave/S-wave velocity at source/receiver (km/s)
        - rho[s|r] = Layer density at source/receiver (kg/m^3)
        - i_0 = Incidence angle at receiver (rad)
        - i_h = Take-off angle at source (rad)
        - dihddel = change in incident angle with distance (rad/rad)
        - p = Horiz. slowness (ray parameter) in layered Earth (not spherical)
    Output parameters:
    """

    # Set constants
    mod_name = 'iasp91'
    del_inc = 0.1  # Increment for del to calculate dihddel
    rad_E = 6371.0  # Earth radius
    # End of constants

    # Make array of phases containing down-going and up-going legs
    # (see TauP manual)
    phases = [phase[0], phase[0].lower()]

    # 1. Get velocities at receiver and source depth ########################
    # Find velocity at source depth. Read in table, then find nearest vel.
    # Find receiver velocity, too.

    # ObsPy path that contains velocity model in tvel format
    vel_mod = ("{:}/{:}.tvel".format(
        os.path.join(
            os.path.dirname(os.path.abspath(obspy.taup.__file__)), "data"),
        mod_name))

    # Import table and find velocities
    vel_table = [
        (float(lyr.split()[0]), float(lyr.split()[1]), float(lyr.split()[2]))
        for n, lyr in enumerate(open(vel_mod)) if n > 1]
    vp_s = [layer for n, layer in enumerate(vel_table) if layer[0] <= z_src
            and vel_table[n+1][0] > z_src and n < len(vel_table)-1][0][1]
    vs_s = [layer for n, layer in enumerate(vel_table) if layer[0] <= z_src
            and vel_table[n+1][0] > z_src and n < len(vel_table)-1][0][2]
    vp_r = vel_table[0][1]
    vs_r = vel_table[0][2]

    # 2. Compute density at source and receiver #############################
    # Coefficients from Nafe-Drake curve (Brocher, 2005 BSSA) for Vp < 8.5 km/s
    # Use linear Birch law approximation for Vp > 8.5 km/s
    if vp_s < 8.5:
        rho_s = ((1.6612 * vp_s - 0.4721 * vp_s ** 2 + 0.0671 * vp_s ** 3
                  - 0.0043 * vp_s ** 4 + 0.00011 * vp_s ** 5) * 1000.0)
    else:
        rho_s = (0.61 + vp_s * 0.328) * 1000.0   # Check this.

    if vp_r < 8.5:
        rho_r = ((1.6612 * vp_r - 0.4721 * vp_r ** 2 + 0.0671 * vp_r ** 3
                  - 0.0043 * vp_r ** 4 + 0.00011 * vp_r ** 5) * 1000.0)
    else:
        rho_r = (0.61 + vp_r * 0.328) * 1000.0   # Check this.

    # 2. Now do free spreading correction ###################################
    # Get first arrival from TauP for incidence angle and ray param
    taup_model = obspy.taup.TauPyModel(model=mod_name)
    arrival = taup_model.get_travel_times(source_depth_in_km=z_src,
                                          distance_in_degree=dist,
                                          phase_list=phases)[0]

    # Estimate dih/dΔ (dihdel) - change of takeoff angle with distance
    if dist - del_inc >= 0:
        dist_1 = dist - del_inc
    else:
        dist_1 = 0.0
    ang_toff_1 = taup_model.get_travel_times(
        source_depth_in_km=z_src, distance_in_degree=dist_1,
        phase_list=phases)[0].takeoff_angle
    dist_2 = dist + del_inc
    ang_toff_2 = taup_model.get_travel_times(
        source_depth_in_km=z_src, distance_in_degree=dist_2,
        phase_list=phases)[0].takeoff_angle
    dihddel = np.abs(np.deg2rad(ang_toff_2 - ang_toff_1) /
                     np.deg2rad(dist_2 - dist_1))

    # Compute angles and slowness. i_0 = incidence angle; i_h = take-off angle
    i_0 = np.deg2rad(arrival.incident_angle)
    i_h = np.deg2rad(arrival.takeoff_angle)

    # Select velocity to use
    if phase == "P":
        v_source = vp_s
        v_rec = vp_r
    elif phase == "S":
        v_source = vs_s
        v_rec = vs_r
    
    # Calculate geometric spreading correction
    g = np.sqrt(
        (rho_s * v_source * np.sin(i_h) * dihddel) /
        (rho_r * v_rec * np.sin(np.deg2rad(dist)) * np.cos(i_0)))

    # 3. Now do free surface correction ######################################
    # Aki & Richards; Kennett (1991, GJI); Yoshimoto et al (1997, PEPI)
    if phase == "P":
        # Convert ray parameter spherical (r.sini / v) -> flat (sini / v)
        p = arrival.ray_param / rad_E  # Spherical -> Flat-Earth ray-param
        q_p0 = np.sqrt(1 / vp_r**2 - p**2)
        q_s0 = np.sqrt(1 / vs_r**2 - p**2)
        # Calculate the denominator for the slowness-dependent quantities
        # Also known as the Rayleigh function (Cerveny)
        gamma = (vp_r**(-2) - 2 * p**2)**2 + (4 * p**2 * q_p0 * q_s0)

        U = np.abs(-2 * vp_r * vs_r**(-2) * q_p0 * (vs_r**(-2) - 2 * p**2)
                   / gamma)
#    elif phase == "SV":
#        U = -2 * vs_r**(-1) * q_s0 * (vs_r**(-2) - 2 * p**2) / gamma
    elif phase == "S":  # Use SH as we are using transverse component
        U = 2.0

    # 4. Now compute the RMS (spherical average) Radiation pattern
    if phase == "P":
        radp = np.sqrt(4.0 / 15.0)
    elif phase == "S":
        radp = np.sqrt(2.0 / 5.0)

    # 5. Now combine into single correction (units: m)
    amom = (g * radp * U) / (4 * np.pi * rho_s * (rad_E * 1000) *
                             (v_source*1000)**3)
    
    return amom