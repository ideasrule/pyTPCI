import subprocess
import numpy as np

pluto_radii = np.load("pluto_radii.npy")

STELLAR_SPEC = "spectra.ini"
a = 2.3e11
k_B = 1.38e-16
Rp = 2.9 * 6.378136e8
m_H = 1.67e-24
hydrogen_frac = 0.92
init_temp = 740
init_mu = 1.23

unit_density = 1e-10
unit_length = Rp
unit_velocity = 1e5
unit_pressure = unit_density * unit_velocity**2

iter_increase_at_step = 50
dt = 0.1

def run_cloudy(global_ind):
    pluto_filename = "data." + str(global_ind).zfill(4) + ".tab"            
    radii, rho, v, P = np.loadtxt(pluto_filename, usecols=(0,2,3,6), unpack=True)
    radii *= unit_length
    rho *= unit_density
    v *= unit_velocity
    #CLOUDY doesn't like velocities that are exactly 0, so set to 1 cm/s
    v[v==0] = 1
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
    script += 'init "spectra.ini"\n'
    script += 'radius {:.6e} linear\n'.format(a)
    script += 'stop depth {:.6e} linear\n'.format(0.9995 * depths.max())
    
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
        script += 'iterate 2\n'
        
    script += 'set dynamics advection length fraction 0.01\n'
    script += 'no molecules\n'
    script += 'element limit off -2\n'
    script += 'metals 0 linear\n'
    script += 'stop temperature linear 5 K\n'
    script += 'turbulence 1 km/sec no pressure\n'
    script += 'double optical depth\n'
    script += 'abundances GASS10 no grains\n'
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
    subprocess.run(["../cloudy/source/cloudy.exe", "-r", prefix], check=False)


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
            f.write("{} {}\n".format(pluto_radii[i], heating_interp[i]))

    

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

for i in range(100):
    t = i * dt
    #t = 10 + (i-100) * dt
    write_heating_file(i)
    run_pluto(i, t)
    run_cloudy(i)
