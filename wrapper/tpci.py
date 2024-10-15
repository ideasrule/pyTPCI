import subprocess
import numpy as np
import sys
import scipy.interpolate
import pdb
import os

# Physical constants
AU = 1.5e13
k_B = 1.38e-16
R_earth = 6.378136e8
R_jup = 11.2 * R_earth
M_earth = 5.97e27
M_jup = 1.898e30
M_sun = 2e33
m_H = 1.67e-24

pluto_radii = np.load("pluto_radii.npy")

# System parameters ---
STELLAR_SPEC = "spectra_hd209458_sf.ini"
a = 0.04723 * AU
Rp = 1.359 * R_jup
Mp = 0.69 * M_jup
M_star = 1.131 * M_sun
gamma = 5./3
tau_frac = 0.1

init_temp = 1320
metallicity = 1
# --------------------

min_temperature = 28 #cloudy limit

if metallicity == 1 or metallicity == 0:
    init_mu = 1.28
    hydrogen_frac = 0.909
elif metallicity == 10:
    init_mu = 1.39
    hydrogen_frac = 0.91
elif metallicity == 30:
    init_mu = 1.68
    hydrogen_frac = 0.90
elif metallicity == 100:
    init_mu = 2.62
    hydrogen_frac = 0.84
else:
    assert(False)

# No need to change
unit_density = 1e-10
unit_length = Rp
unit_velocity = 1e5
unit_pressure = unit_density * unit_velocity**2
unit_time = unit_length / unit_velocity

advec_turnon = 100


def write_params_header(filename="params.h"):
    ''' Write PLUTO params file'''
    f = open(filename, "w")
    f.write("#define  PLANET_RADIUS   ({} * CONST_Rearth)\n".format(Rp/R_earth))
    f.write("#define  PLANET_MASS     ({} * CONST_Mearth)\n".format(Mp/M_earth))
    f.write("#define  MU              {}\n".format(init_mu))
    f.write("#define  INIT_TEMP       {}\n".format(init_temp))
    f.write("#define  STAR_MASS       ({} * CONST_Msun)\n".format(M_star/M_sun))
    f.write("#define  SEMIMAJOR_AXIS  ({} * CONST_au)\n".format(a/AU))
    f.write("#define  n_SURF          1e14\n")
    f.write("#define  P_SURF          (n_SURF * CONST_kB * INIT_TEMP)\n")
    f.write("#define  RHO_SURF        (P_SURF * (MU * CONST_mH) / (CONST_kB * INIT_TEMP))\n")
    f.write("#define  BETA            (CONST_G * PLANET_MASS * (MU * CONST_mH) / (PLANET_RADIUS * CONST_kB * INIT_TEMP))\n")
    f.write("#define  UNIT_DENSITY    {}\n".format(unit_density))
    f.write("#define  UNIT_LENGTH     {:.3e}\n".format(Rp))
    f.write("#define  UNIT_VELOCITY   {:.3e}\n".format(unit_velocity))
    f.close()

    
def run_cloudy(global_ind, t, last_cloudy_iter):
    ''' RUN CLOUDY from previous PLUTO data'''
    pluto_filename = "data." + str(global_ind).zfill(4) + ".tab"        
    radii, rho, v, P = np.loadtxt(pluto_filename, usecols=(0,2,3,6), unpack=True)
    radii *= unit_length
    rho *= unit_density
    v *= unit_velocity
    P *= unit_pressure
    depths = np.max(radii) - radii

    if last_cloudy_iter > 0:
        prev_cloudy_file = "cl_data." + str(last_cloudy_iter).zfill(4) + ".over.tab"
        cloudy_depths, cloudy_mu = np.loadtxt(prev_cloudy_file, usecols=(0,3), unpack=True)
        mu = np.interp(depths, cloudy_depths, cloudy_mu)
    else:
        mu = np.ones(len(depths)) * init_mu
        
    Ts = (P * mu * m_H) / (k_B * rho)
    script = ''
    script += 'CMB\n'
    script += 'cosmic rays background\n'
    script += 'init "{}"\n'.format(STELLAR_SPEC)
    script += 'stop depth {:.6e} linear\n'.format(depths.max())
    script += 'illumination angle 66 deg\n'

    n_H = rho / (init_mu * m_H) * hydrogen_frac # density
    script += 'dlaw table depth linear\n'
    for i in range(len(depths)-1, -1, -1):
        script += '{:.6e} {:.3e}\n'.format(depths[i], n_H[i])
    # Fake point to extend range of table
    script += '{:.6e} {:.3e}\n'.format(1.01*depths[0], n_H[0])
    script += 'end of dlaw\n'
    
    script += 'tlaw table depth linear\n' # temperature
    for i in range(len(depths)-1, -1, -1):
        if Ts[i] < min_temperature: # fix to improve CLOUDY stability
            print(f"WARNING: temperature {Ts[i]} below min {min_temperature}")
            Ts[i] = min_temperature
            
        script += '{:.6e} {}\n'.format(depths[i], Ts[i])    
    script += '{:.6e} {:.3e}\n'.format(1.01*depths[0], Ts[0])
    script += 'end of tlaw\n'

    if t > advec_turnon:
        script += 'wind advection table depth linear\n' # wind velocities
        for i in range(len(depths)-1, -1, -1):
            # Positive (inward) velocities crash CLOUDY
            script += '{:.6e} {:.3e}\n'.format(depths[i], min(-1e-10, -v[i]))
        script += '{:.6e} {:.3e}\n'.format(1.01*depths[0], -1e-10)
        script += 'end of velocity table\n'
        script += 'set dynamics advection length fraction 0.01\n'
        #script += 'set dynamics advection length 9\n'
        script += 'iterate to convergence max=40\n'
    else:
        script += 'iterate 2\n'

    #script += 'failures 100\n'
    #script += 'set nend 2800\n'
    script += 'no molecules\n'
    script += 'stop temperature linear 5 K\n'
    script += 'turbulence 1 km/sec no pressure\n'
    script += 'double optical depth\n' # simulate planet's surface
    script += 'element limit off -5\n'
    script += 'abundances GASS10 no grains\n'
    script += 'metals {} linear\n'.format(metallicity)
    script += 'print short\n'
    script += 'print line faint -2 log\n'
    prefix = "cl_data." + str(global_ind).zfill(4)
    script += 'set save prefix "{}"\n'.format(prefix + ".")
    script += 'set save hash ""\n'
    script += 'save overview "all_iters_over.tab"\n'
    script += 'save overview "over.tab" last\n'
    script += 'save pressure "pres.tab" last\n'
    script += 'save wind "wind.tab"\n'
    #script += 'save dynamics advection "dyna.tab" last\n'
    script += 'save continuum "continuum.tab" last units Angstrom\n'
    script += 'save cooling "cool.tab" last\n'
    script += 'save ages "ages.tab" last\n'
    script += 'save species labels all "labels.tab"\n'
    script += 'save species populations "pops.tab" all\n'
    script += 'save species populations "He.tab" "He[1:2]"\n'
    script += 'save species populations "H.tab" "H[1:2]"\n'
    script += 'save species energies "energies.tab" all\n'

    with open(prefix + ".in", "w") as f:
        f.write(script)

    print("Running cloudy", global_ind)
    result = subprocess.run(["../cloudy/source/cloudy.exe", "-r", prefix], check=False)
    if result.returncode != 0:
        print("Warning: something went wrong in CLOUDY.  You can ignore this in the early timesteps")
        

def write_heating_file(global_ind, output_file="curr_heat.dat"):
    ''' Write heating file curr_heat.dat from CLOUDY to transfer into PLUTO'''
    cloudy_file = "cl_data." + str(global_ind).zfill(4) + ".over.tab"
    data = np.loadtxt(cloudy_file, usecols=(0,1,5,6))
    data = data[::-1]
    depth, rho, heating, cooling = data.T
    net_heating = heating - cooling
    net_heating /= rho
    cloudy_radii = np.max(pluto_radii) - depth / unit_length

    interpolator = scipy.interpolate.interp1d(cloudy_radii, net_heating, bounds_error=False, fill_value=0)
    pluto_widths = np.gradient(pluto_radii)
    heating_interp = interpolator(pluto_radii)
    
    for i in range(len(pluto_radii)):
        cond = (cloudy_radii >= pluto_radii[i] - pluto_widths[i] / 2) & (cloudy_radii < pluto_radii[i] + pluto_widths[i] / 2)
        #print(i, np.sum(cond))
        if np.sum(cond) > 1:
            heating_interp[i] = np.mean(net_heating[cond])
        

    with open(output_file, "w") as f:
        for i in range(len(pluto_radii)):
            f.write("{} {:.6e}\n".format(pluto_radii[i], heating_interp[i]))
    

def run_pluto(global_ind, t, template_file="pluto_template.ini", ini_file="pluto.ini"):
    ''' Run PLUTO '''
    with open(template_file, "r") as f:
        template = f.read()
    config = template.replace("TSTOP", str(t+dt))
    with open(ini_file, "w") as f:
        f.write(template.replace("TSTOP", str(t+dt)))

    command = ["./pluto"]
    if t > 0:
        command += ["-h5restart", "{}".format(global_ind)]
    print("Running pluto command: ", command)
    f = open("pluto_log.txt", "a")
    subprocess.run(command, check=True, stdout=f)
    f.close()

    
def get_max_change(global_ind, last_cloudy_iter):
    ''' Find maximum relative difference of density and pressure '''
    if global_ind == 0:
        return False
    
    prev_pluto_filename = "data." + str(last_cloudy_iter).zfill(4) + ".tab"
    radii, rho, P = np.loadtxt(prev_pluto_filename, usecols=(0,2,6), unpack=True)
    curr_pluto_filename = "data." + str(global_ind).zfill(4) + ".tab"
    radii, rho2, P2 = np.loadtxt(curr_pluto_filename, usecols=(0,2,6), unpack=True)

    max_rho_rel_diff = np.max(np.abs(rho2 / rho - 1))
    max_P_rel_diff = np.max(np.abs(P2 / P - 1))
    max_rel_diff = max(max_rho_rel_diff, max_P_rel_diff)
    print("Max rel change: ", max_rel_diff)
    return max_rel_diff

def delete_pluto_files(last_cloudy_iter, curr_cloudy_iter):
    '''Delete extraneous PLUTO files with no corresponding CLOUDY runs'''
    for i in range(last_cloudy_iter+1, curr_cloudy_iter):
        prefix = "data." + str(i).zfill(4)
        try:
            os.remove(f"{prefix}.dbl")
            os.remove(f"{prefix}.dbl.xmf")
            os.remove(f"{prefix}.dbl.h5")
            os.remove(f"{prefix}.tab")
        except:
            pass

write_params_header()
#subprocess.run(["make", "clean"], check=True)
subprocess.run(["make"], check=True)

# [start, t, dt]
if len(sys.argv) == 4:
    global_ind = int(sys.argv[1])
    t = float(sys.argv[2])
    dt = float(sys.argv[3])
    last_cloudy_iter = global_ind
elif len(sys.argv) == 1:
    global_ind = 0
    t = 0
    dt = 1e-2
    last_cloudy_iter = 0
else:
    print("Give no arguments to start from beginning.  Give timestep and time to resume.")
    assert(False)

max_t = 1000 # run pyTPCI until max_t
log_f = open("tpci_log.txt", "a")

# main loop
while t < max_t:
    log_f.write("Starting loop for step {}, t={}, dt={}\n".format(
        global_ind, round(t,5), round(dt,5)))
    log_f.flush()
    run_pluto(global_ind, t)
    max_rho_change = get_max_change(global_ind, last_cloudy_iter)
    if max_rho_change > 0.1:
        #delete_pluto_files(last_cloudy_iter, global_ind)
        run_cloudy(global_ind, t, last_cloudy_iter)
        write_heating_file(global_ind)
        last_cloudy_iter = global_ind
    
    t += dt
    global_ind += 1

log_f.close()
