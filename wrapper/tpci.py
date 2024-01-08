import subprocess
import numpy as np
import sys

#Physical constants
AU = 1.5e13
k_B = 1.38e-16
R_earth = 6.378136e8
M_earth = 5.97e27
m_H = 1.67e-24

pluto_radii = np.load("pluto_radii.npy")

#System parameters
STELLAR_SPEC = "spectra.ini"
a = 0.06 * AU
Rp = 2.9 * R_earth
Mp = 11 * M_earth
M_sun = 2e33
M_star = 0.75 * M_sun

init_temp = 740
metallicity = 0

if metallicity == 1 or metallicity == 0:
    init_mu = 1.23
    hydrogen_frac = 0.92
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

#Probably no need to change
unit_density = 1e-10
unit_length = Rp
unit_velocity = 1e5
unit_pressure = unit_density * unit_velocity**2

iter_increase_at_step = 50


def write_params_header(filename="params.h"):
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
        
def run_cloudy(global_ind):
    pluto_filename = "data." + str(global_ind).zfill(4) + ".tab"            
    radii, rho, v, P = np.loadtxt(pluto_filename, usecols=(0,2,3,6), unpack=True)
    radii *= unit_length
    rho *= unit_density
    v *= unit_velocity
    P *= unit_pressure
    depths = np.max(radii) - radii

    if global_ind > 0:
        prev_cloudy_file = "cl_data." + str(global_ind-1).zfill(4) + ".over.tab"
        cloudy_depths, cloudy_mu = np.loadtxt(prev_cloudy_file, usecols=(0,3), unpack=True)
        mu = np.interp(depths, cloudy_depths, cloudy_mu)
    else:
        mu = np.ones(len(depths)) * init_mu
        
    Ts = (P * mu * m_H) / (k_B * rho)
    script = ''
    script += 'CMB\n'
    script += 'cosmic rays background\n'
    script += 'init "{}"\n'.format(STELLAR_SPEC)
    script += 'radius {:.6e} linear\n'.format(a)
    script += 'stop depth {:.6e} linear\n'.format(0.9995 * depths.max())
    script += 'illumination angle 66 deg\n'
    
    script += 'dlaw table depth linear\n'
    for i in range(len(depths)-1, -1, -1):
        script += '{:.6e} {:.3e}\n'.format(depths[i], rho[i] / (init_mu * m_H) * hydrogen_frac)
    script += 'end of dlaw\n'
    
    script += 'tlaw table depth linear\n'
    for i in range(len(depths)-1, -1, -1):
        script += '{:.6e} {}\n'.format(depths[i], Ts[i])
    script += 'end of tlaw\n'

    script += 'wind advection table depth linear\n'
    for i in range(len(depths)-1, -1, -1):
        #Positive (inward) velocities crash CLOUDY for some reason
        script += '{:.6e} {:.3e}\n'.format(depths[i], min(-1, -v[i]))
    script += 'end of velocity table\n'

    if global_ind > iter_increase_at_step:
        script += 'iterate 20\n'
    else:
        script += 'iterate 1\n'
        
    script += 'set dynamics advection length fraction 0.01\n'
    script += 'no molecules\n'
    script += 'stop temperature linear 5 K\n'
    script += 'turbulence 1 km/sec no pressure\n'
    script += 'double optical depth\n'
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
    script += 'save dynamics advection "dyna.tab" last\n'
    script += 'save continuum "continuum.tab" last units Angstrom\n'
    script += 'save cooling "cool.tab" last\n'
    script += 'save ages "ages.tab" last\n'
    script += 'save species labels all "labels.tab"\n'
    script += 'save species populations "pops.tab" all\n'
    script += 'save species populations "He.tab" "He[1:2]"\n'
    script += 'save species energies "energies.tab" all\n'

    with open(prefix + ".in", "w") as f:
        f.write(script)

    print("Running cloudy", global_ind)
    result = subprocess.run(["../cloudy/source/cloudy.exe", "-r", prefix], check=False)
    if result.returncode != 0:
        print("Warning: something went wrong in CLOUDY.  You can ignore this in the early timesteps")


def write_heating_file(global_ind, output_file="curr_heat.dat"):
    if global_ind == 0:
        heating_interp = np.zeros(len(pluto_radii))
    else:
        cloudy_file = "cl_data." + str(global_ind-1).zfill(4) + ".over.tab"
        depth, rho, heating = np.loadtxt(cloudy_file, usecols=(0,1,5), unpack=True)
        heating /= rho
        cloudy_radii = 1 + (np.max(depth) - depth) / unit_length
        heating_interp = np.interp(pluto_radii, cloudy_radii[::-1], heating[::-1])

    with open(output_file, "w") as f:
        for i in range(len(pluto_radii)):
            f.write("{} {:.6e}\n".format(pluto_radii[i], heating_interp[i]))

    

def run_pluto(global_ind, t, template_file="pluto_template.ini", ini_file="pluto.ini"):
    with open(template_file, "r") as f:
        template = f.read()
    config = template.replace("TSTOP", str(t+dt))
    with open(ini_file, "w") as f:
        f.write(template.replace("TSTOP", str(t + dt)))

    command = ["./pluto"]
    if t > 0:
        command += ["-restart", "{}".format(global_ind)]
    print("Running pluto command: ", command)
    subprocess.run(command, check=True)

def is_converged(global_ind, max_frac_diff=0.02):
    if global_ind == 0:
        return False
    
    prev_pluto_filename = "data." + str(global_ind-1).zfill(4) + ".tab"
    radii, rho = np.loadtxt(prev_pluto_filename, usecols=(0,2), unpack=True)
    curr_pluto_filename = "data." + str(global_ind).zfill(4) + ".tab"
    radii, rho2 = np.loadtxt(curr_pluto_filename, usecols=(0,2), unpack=True)

    max_rho_rel_diff = np.max(np.abs(rho2 / rho - 1))
    print("Max rel change in rho: ", max_rho_rel_diff)
    return max_rho_rel_diff < max_frac_diff
    

write_params_header()
#subprocess.run(["make", "clean"], check=True)
subprocess.run(["make"], check=True)

if len(sys.argv) == 3:
    global_ind = int(sys.argv[1])
    t = float(sys.argv[2])
elif len(sys.argv) == 1:
    global_ind = 0
    t = 0
else:
    print("Give no arguments to start from beginning.  Give timestep and time to resume.")
    assert(False)

dt = 0.1
max_t = 100

log_f = open("tpci_log.txt", "w")

while t < max_t:
    log_f.write("Starting loop for step {}, t={}, dt={}\n".format(
        global_ind, round(t,5), round(dt,5)))
    log_f.flush()
    write_heating_file(global_ind)
    run_pluto(global_ind, t)
    run_cloudy(global_ind)
    
    if is_converged(global_ind):
        dt *= 2
        print("Doubling dt to {}".format(dt))

    t += dt
    global_ind += 1

log_f.close()
