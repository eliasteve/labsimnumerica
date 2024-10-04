import os
import subprocess

def run_simulation(
    folder, #Location of the folder NSL_SMIULATOR
    # Input properties as they will be written in input.dat
    sim_type=0, 
    J = 1.0,
    H = 0.0,
    restart=0,
    temperature=1.1,
    n_part=108,
    density=0.8,
    r_cut=2.5,
    delta=0.0005,
    n_blocks=10,
    n_steps=1000,
    n_eqsteps=0,
    tail=0
):
  """
  Writes input.dat and then runs a simulation with the specified
  parameters.

  Args:
    folder: path of the NSL_SIMULATOR folder
    all other args: arguments as they will be written in input.dat
  """
  if (sim_type in {2, 3}):
      input_as_string = f"SIMULATION_TYPE        {sim_type} {J} {H}\n"
  else:
      input_as_string = f"SIMULATION_TYPE        {sim_type}\n"
  input_as_string += f"RESTART                {restart}\n"
  #Putting in all the digits of the temperature causes a bug in the sim for
  #Ising.
  input_as_string += f"TEMP                   {temperature:1.4f}\n"
  input_as_string += f"NPART                  {n_part}\n"
  input_as_string += f"RHO                    {density}\n"
  input_as_string += f"R_CUT                  {r_cut}\n"
  input_as_string += f"DELTA                  {delta}\n"
  input_as_string += f"NBLOCKS                {n_blocks}\n"
  input_as_string += f"NSTEPS                 {n_steps}\n"
  input_as_string += f"NEQSTEPS               {n_eqsteps}\n"
  input_as_string += f"TAIL                   {tail}\n"
  input_as_string += "\n"
  input_as_string += "ENDINPUT\n"

  with open(folder + "/INPUT/input.dat", "w") as f:
      f.write(input_as_string)
  cwd = os.getcwd()
  try:
    os.chdir(folder + "/SOURCE")
    res = subprocess.run(["./simulator.exe"])
    os.chdir(cwd)
  #To avoid interrupting the program during the simulation and not
  #changing the current directory to the previous value
  except KeyboardInterrupt:
    os.chdir(cwd)

def set_properties(folder, properties, gofr_points = None):
  """
  Writes properties.dat with the properties to be measured in the
  simulation.

  Args:
    folder: path of the NSL_SIMULATOR folder
    properties: list of the properties to be measured. Do not need to be
      capitalized.
    gofr_points: how many points to use in the evaluation of the radial
      distribution function
  """
  with open(folder + "/INPUT/properties.dat", "w") as f:
    for prop in properties:
      if prop.upper() == "GOFR" and gofr_points != None:
        f.write(prop.upper() + f" {gofr_points}\n")
      else:
        f.write(prop.upper() + "\n")
    f.write("\nENDPROPERTIES\n")

def set_config_file(folder, mode):
  """
  Writes the appropriate config file for the starting configuration (either
  for Lennard-Jones of Ising system)

  Args:
    folder: path of the NSL_SIMULATOR folder
    mode: one of "ising" (for writing the config appropriate for Ising
      systems) or "lj" (for writing the config appropriate for 
      Lennard-Jones systems)
  """
  if mode == "ising":
    subprocess.run(["cp", folder+"/INPUT/CONFIG/config.ising", folder+"/INPUT/CONFIG/config.xyz"])
  elif mode == "lj":
    subprocess.run(["cp", folder+"/INPUT/CONFIG/config.fcc", folder+"/INPUT/CONFIG/config.xyz"])
  else:
    raise ValueError("Unrecognized mode (options are ising or lj)")

def remove_output(folder):
  """
  Removes the simulation output from the NSL_SIMULATOR folder. Analogous to
  running "make remove".

  Args:
    folder: path of the NSL_SIMULATOR folder
  """ 
  subprocess.run([f"rm {folder}/OUTPUT/*.*"], shell=True) #shell=True for filename wildcard
  subprocess.run([f"rm {folder}/OUTPUT/CONFIG/*.*"], shell=True)
