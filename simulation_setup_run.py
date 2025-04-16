import re
import sys
from os import system
import os
from datetime import datetime
import shutil
import socket

try:
    actions = sys.argv
    for action in actions[2::]:
        if action not in ["compile", "make_clean", "setup", "run", "clean_output", "clean_setup_log"]:
            raise ValueError

except:
    print(
        "============================================================================================================"
    )
    print(
        "Usage:\npython3 simulation_setup_run.py control_filepath compile setup run clean_output clean_setup_log (at least one of last five)"
    )
    print(
        "\ncontrol_filepath is the location of the control file relative to simulation_setup_run.py"
    )
    print("\n==========================================================")
    print("Actions: ")
    print("==========================================================")
    print("\ncompile: compile src code and create executable file")
    print("\nmake_clean: remove all .mod and .o files")
    print(
        "\nsetup: read control.dat and generates ind3dmg.dat, ind3dmg_LSM.dat, nproc.dat and queue file"
    )
    print("\nrun: submit the job to be executed")
    print("\nclean_output: remove all files in output folder")
    print("\nclean_setup_log: remove all files in setup_log folder")

    print(
        "\n \nExample:\n \npython3 simulation_setup_run.py ../control.dat compile run"
    )
    print(
        "============================================================================================================"
    )

    sys.exit(1)

# Date and time for log_files
now = datetime.now()
dt_string = now.strftime("%Y_%m_%d_%H_%M_%S")
current_working_directory = os.getcwd()


def execute(command):
    print("> " + command)
    system(command)

def make_clean(parameters):
    p = parameters
    src_directory = f"{current_working_directory}/{p['source_code']}"
    os.chdir(src_directory)
    execute(f"make clean")
    os.chdir(current_working_directory)

def compile(parameters):
    p = parameters

    # Get the hostname of the server
    server_name = socket.gethostname()
   
    src_directory = f"{current_working_directory}/{p['source_code']}"
    os.chdir(src_directory)
    #execute(f"make")

    if server_name == 'Dellemc':
        execute(f"make fondecyt")
    elif server_name == 'leftraru2':
        execute(f"make leftraru")
    else:
        # Default case if the server name is not recognized
        print(f"Unknown server: {server_name}. The compilation will be carried out using make only")
        execute(f"make")

    execute(f"mv channel ../")
    os.chdir(current_working_directory)


def get_data():

    # This routine reads every line of the control file. If in some line there is a "=" symbol, it adds an entry to the dictionary
    # "parameters" with the name at the left-hand side of the equal symbol with the value of the right-hand side.
    # It returns the dictionary "parameters"
    control_file_path = actions[1]
    # empty dictionary to save the parameters
    parameters = {}

    # open the control file

    print(f"\nReading control.dat ...")

    f = open(control_file_path, "r")

    # loop reading the entries of the control file
    for l in f:
        parameter_entry = l.split("=")
        if len(parameter_entry) > 1:  # it means there is a '=' in
            # line and therefore it is an entry and not a header
            parameters[parameter_entry[0].strip()] = parameter_entry[1].strip()
            # print(parameter_entry[0].strip(), " = ", parameter_entry[1].strip())
    f.close()

    print(f"\ncontrol.dat reading completed!")

    return parameters


def write_ind3dmg(parameters):

    # to shorten the coding lines
    p = parameters 

    # time_steps_to_save is a list with all the time steps to be saved for debugging
    time_steps_to_save_string     = f"{p['NS_time_steps_to_save']}" # string: "t1, t2 ... tn"
    
    if(time_steps_to_save_string != "0"): # there are time steps to monitor   

        # Now I need to split the string by commas and strip it by spaces
        time_steps_to_save_list       = time_steps_to_save_string.split(',')
        # convert it to a list. Example: time_steps_to_save_debugging = [1, 100, 1000, 2000] 
        time_steps_to_save_debugging  = [ int(tsteps.strip(' ')) for tsteps in time_steps_to_save_list ]
        ntime_steps_to_save_debugging = len(time_steps_to_save_debugging)

    else: # there are no time steps to monitor
        ntime_steps_to_save_debugging = 0

    bcs = {
        "interface": "0",
        "wall": "1",
        "symmetry_plane": "2",
        "freestream": "3",
        "inflow": "4",
        "MOC": "5",
        "periodic": "6",
        "outflow": "8",
        "slip_wall": "9",
    }

    # If ind3dmg.dat exists, I rename it with the date and time info
    if os.path.exists("ind3dmg.dat"):
        os.rename("ind3dmg.dat", f"setup_log/ind3dmg_old_{dt_string}.dat")

    print(f"\nGenerating ind3dmg.dat ...")

    f = open("ind3dmg.dat", "a")
    f.write(f"1 1\n")
    f.write(f"{p['imax']} {p['jmax']} {p['kmax']}\n")
    f.write(f"0 0 0\n")
    f.write(
        f"{p['total_time_steps']} {p['dt']} {p['pseudo_time_min_iterations']} -2.0 -7.5\n"
    )
    f.write(f"1\n")
    f.write(f"{p['pseudo_time_max_iterations']}\n")
    f.write(f"{p['Reynolds_number']} {p['beta']} 1 {p['usolu_save_iterations']}\n")
    f.write(f"{p['CFL_number']} {p['Von_Neumann_number']}\n")
    f.write(f"4 0.25 0.33333333333 0.5 1.0\n")
    f.write(f"0.065 0.065 0.065 0.065\n")
    f.write(f"2 0.65 1.0E-08 1.0E-06\n")
    f.write(f"{p['CFL_turbulence_model']} {p['Von_Neumann_turbulence_model']}\n")
    f.write(f"{p['eps_x']} {p['eps_y']} {p['eps_z']}\n")
    f.write(f"0.2 0.2 0.2\n")
    f.write(f"{p['solu_restart_update_iterations']}\n")
    f.write(f"{p['pressure_dissipation_coefficient']}\n")
    f.write(f"{p['monitoring_points']}\n")
    f.write(
        f"{p['point1_coordinates'].split(',')[0]} {p['point1_coordinates'].split(',')[1]} {p['point1_coordinates'].split(',')[2]} {p['point1_coordinates'].split(',')[3]}\n"
    )
    f.write(
        f"{bcs[p['i1']]} {bcs[p['im']]} {bcs[p['j1']]} {bcs[p['jm']]} {bcs[p['k1']]} {bcs[p['km']]}\n"
    )
    f.write(f"1 1 2 2 3 3\n")
    f.write(f"1 1 1 1\n")
    f.write(f"{p['call_solver_daf']}\n")
    f.write(f"{p['hydraulic_mode']}\n")
    f.write(f"{p['dynamic_dtau']}\n")
    f.write(f"{p['non_slip_wall_blanking']}\n")
    f.write(f"{ntime_steps_to_save_debugging}\n")

    if( ntime_steps_to_save_debugging >0 ):
        for n in range(ntime_steps_to_save_debugging):
            f.write(f"{time_steps_to_save_debugging[n]}\n")    
    else:
        f.write(f"0\n")

    #if(f"{p['blanking_zones']}" != "0"):
    #    f.write(f"{p['blanking_zones']}\n")
    #
    #    for n in range(int(p['blanking_zones'])):
    #        f.write(f"0 0 0 0 0 0\n")
    #        f.write(f"{p['iblk1']} {p['iblk2']} \n")
    #        f.write(f"{p['jblk1']} {p['jblk2']} \n")
    #        f.write(f"{p['kblk1']} {p['kblk2']} \n")


    if f"{p['blanking_zones']}" != "0":

        f.write(f"{p['blanking_zones']}\n")

        # Remove the brackets
        iblk1_str_clean = p['iblk1'].strip("[]").strip()
        iblk2_str_clean = p['iblk2'].strip("[]").strip()
        jblk1_str_clean = p['jblk1'].strip("[]").strip()
        jblk2_str_clean = p['jblk2'].strip("[]").strip()
        kblk1_str_clean = p['kblk1'].strip("[]").strip()
        kblk2_str_clean = p['kblk2'].strip("[]").strip()
                    
        iblk1_ints_array = [int(num.strip()) for num in iblk1_str_clean.split(",")]
        iblk2_ints_array = [int(num.strip()) for num in iblk2_str_clean.split(",")]
        jblk1_ints_array = [int(num.strip()) for num in jblk1_str_clean.split(",")]
        jblk2_ints_array = [int(num.strip()) for num in jblk2_str_clean.split(",")]
        kblk1_ints_array = [int(num.strip()) for num in kblk1_str_clean.split(",")]
        kblk2_ints_array = [int(num.strip()) for num in kblk2_str_clean.split(",")]

        for n in range(int(p['blanking_zones'])):

            f.write("0 0 0 0 0 0\n")
            f.write(f"{iblk1_ints_array[n]} {iblk2_ints_array[n]}\n")
            f.write(f"{jblk1_ints_array[n]} {jblk2_ints_array[n]}\n")
            f.write(f"{kblk1_ints_array[n]} {kblk2_ints_array[n]}\n")

    f.write(f"0\n")
    f.write(f"0\n")
    f.write(f"0\n")
    f.write(f"0\n")

    f.close()
    print(f"\nind3dmg.dat generated!")


def write_ind3dmg_LSM(parameters):

    p = parameters

    # time_steps_to_save is a list with all the time steps to be saved for debugging
    time_steps_to_save_string     = f"{p['LS_time_steps_to_save']}" # string: "t1, t2 ... tn"
    
    if(time_steps_to_save_string != "0"): # there are time steps to monitor   

        # Now I need to split the string by commas and strip it by spaces
        time_steps_to_save_list       = time_steps_to_save_string.split(',')
        # convert it to a list. Example: time_steps_to_save = [1, 100, 1000, 2000] 
        time_steps_to_save_debugging  = [ int(tsteps.strip(' ')) for tsteps in time_steps_to_save_list ]
        ntime_steps_to_save_debugging = len(time_steps_to_save_debugging)
    
    else: # there are no time steps to monitor
        ntime_steps_to_save_debugging = 0

    bcs = {
        "absorption": "0",
        "extrapolation": "1",
        "advection": "2",
    }
    extrp_opt = {
        "upwind": "1",
        "central_difference": "0",
    }
    narrowband_options = {
        "off": "0",
        "on": "1",
    }
    advection_scheme_options = {
        "UPWIND": "0",
        "WENO3": "1",
    }

    sussman_correction_method_options = {
        "point_constant": "0",
        "simpson_rule": "1",
    }

    # If ind3dmg_LSM.dat exists, I rename it with the date and time info
    if os.path.exists("ind3dmg_LSM.dat"):
        os.rename("ind3dmg_LSM.dat", f"setup_log/ind3dmg_LSM_old_{dt_string}.dat")

    print(f"\nGenerating ind3dmg_LSM.dat ...")
    f = open("ind3dmg_LSM.dat", "a")
    f.write(f"{p['epsilon']} {p['level_set_advection_iterations']}\n")
    f.write(f"0.0\n")
    f.write(
        f"{p['rho_air/rho_ref']} {p['mu_air/mu_ref']} {p['pressure_dissipation_LSM_air']} {p['pressure_dissipation_LSM_ref']}\n"
    )
    f.write(f"{p['Froude_number']} {p['Weber_number']}\n")
    f.write(f"0 700 0.005 0\n")
    f.write(
        f"{p['reinitialization_iterations']} {p['reinitialization_dtau']} {p['reinitilization_iteration_frequency']}\n"
    )
    f.write(f"{p['phi_save_iterations']}\n")
    f.write(f"0\n")
    f.write(
        f"{bcs[p['i1_LSM']]} {bcs[p['im_LSM']]} {bcs[p['j1_LSM']]} {bcs[p['jm_LSM']]} {bcs[p['k1_LSM']]} {bcs[p['km_LSM']]}\n"
    )
    f.write(f"{extrp_opt[p['extrapolation_option']]}\n")
    f.write(f"{narrowband_options[p['narrow_band']]} {p['narrow_band_thickness']}\n")
    f.write(f"{advection_scheme_options[p['advection_scheme']]}\n")
    f.write(f"{p['number_of_obstacles']}\n")
    f.write(f"{sussman_correction_method_options[p['sussman_correction_method']]}\n")
    f.write(f"{p['call_levelsetmethod']}\n")
    f.write(f"{p['call_reinitialisation']}\n")
    f.write(f"{p['hybrid_reinitialisation']}\n")
    f.write(f"{p['sweep_lsqm']}\n")
    f.write(f"{p['radius_lsqm']}\n")
    f.write(f"{p['ConvergenceToleranceGeomReini']}\n")
    f.write(f"{p['TotalVolumeComputation']}\n")
    f.write(f"{p['GlobalMassCorrection']}\n")
    f.write(f"{ntime_steps_to_save_debugging}\n")
    if( ntime_steps_to_save_debugging >0 ):
        for n in range(ntime_steps_to_save_debugging):
            f.write(f"{time_steps_to_save_debugging[n]}\n")    
    f.write(f"{p['epsReinitialisation']}\n")
    f.write(f"{p['OrderLSAdvectionBoundaries']}\n")
    f.write(f"{p['OrderReinitialisationBoundaries']}\n")
    f.write(f"{p['ENOBCReinitialisation']}\n")
    f.write(f"{p['BigPhi']}\n")
    f.write(f"{p['limit_ghost_velocities']}\n")
    f.write(f"{p['zero_pressure_fs']}\n")

    print(f"\nind3dmg_LSM.dat generated!")

    f.close()


def write_nproc(parameters):
    p = parameters

    # If nproc.dat exists, I rename it with the date and time info
    if os.path.exists("nproc.dat"):
        os.rename("nproc.dat", f"setup_log/nproc_old_{dt_string}.dat")

    print(f"\nGenerating nproc.dat ...")
    f = open("nproc.dat", "a")
    f.write(f"{p['number_of_processes']}\n")
    f.write(
        f"{p['i_direction_processes']} {p['j_direction_processes']} {p['k_direction_processes']}"
    )

    f.close()
    print(f"\nnproc.dat generated!")


def write_directory(parameters):
    p = parameters

    # If directory file exists, I rename it with the date and time info
    if os.path.exists("directory"):
        os.rename("directory", f"setup_log/directory_old_{dt_string}")

    current_directory = os.getcwd()
    directory = f"{current_directory}/{p['output_directory']}"


    # if the output_directory folder doesn't exist, we create it

    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Create subfolders
    combined_dir      = f"{directory}combined/"
    debug_dir         = f"{directory}debug/"
    debugfiles_dir    = f"{directory}debugfiles/"
    solufiles_dir     = f"{directory}solufiles/"
    solus_dir         = f"{directory}solus/"
    surface_dir       = f"{directory}surface/"
    usolufiles_dir    = f"{directory}usolufiles/"

    decofiles_dir    = f"{current_directory}/{p['decofiles_directory']}"

    filestat_dir      = f"{directory}filestat"


    if not os.path.exists(combined_dir):
        os.makedirs(combined_dir)

    if not os.path.exists(debug_dir):
        os.makedirs(debug_dir)

    if not os.path.exists(debugfiles_dir):
        os.makedirs(debugfiles_dir)

    if not os.path.exists(solufiles_dir):
        os.makedirs(solufiles_dir)

    if not os.path.exists(solus_dir):
        os.makedirs(solus_dir)

    if not os.path.exists(surface_dir):
        os.makedirs(surface_dir)

    if not os.path.exists(usolufiles_dir):
        os.makedirs(usolufiles_dir)

    if not os.path.exists(decofiles_dir):
        os.makedirs(decofiles_dir)

    if not os.path.exists(filestat_dir):
        print(f"\nGenerating filestat ...")
        f = open("filestat", "w")
        f.write('0.0 \t 0')
        f.close()
        print(f"\noutput directory path file generated!")


    print(f"\nGenerating output directory path file ...")
    f = open("directory", "w")
    f.write(directory)
    f.close()
    print(f"\noutput directory path file generated!")


def write_queue(parameters):
    p = parameters
    control_file_path = actions[1]
    # empty dictionary to save the parameters
    queue_file_lines = []

    # open the control file
    f = open(control_file_path, "r")

    # loop reading the entries of the control file

    reading_queue_parameters = False

    for l in f:
        if "#!" in l or reading_queue_parameters == True:
            reading_queue_parameters = True
            queue_file_lines.append(l)

    f.close()

    # If queue file exists, I rename it with the date and time info
    if os.path.exists(p["queue_filename"]):
        os.rename(
            p["queue_filename"], f"setup_log/{p['queue_filename']}_old_{dt_string}"
        )

    print(f"\nGenerating queue file ...")
    f = open(p["queue_filename"], "a")

    for l in queue_file_lines:
        f.write(l)
    f.close()
    print(f"\nqueue file generated!")


def run(parameters):
    p = parameters

    # Get the hostname of the server
    server_name = socket.gethostname()

    # When I run the simulation, it creates a dated copy of control file into output
    # directory

    #current_directory = os.getcwd()
    directory = f"{current_working_directory}/{p['output_directory']}"
    control_file_path = actions[1]
    shutil.copy(control_file_path, f"{directory}/{control_file_path.split('.dat')[0]}_{dt_string}.dat")


    # generation of simulation files before submitting the job for execution

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # GRID FILE
    # if grid file doesn't exist, I copy the grid file from control.dat as "grid"
    #if not os.path.exists(f"{current_working_directory}/grid"):

    shutil.copy(f"{current_working_directory}/{p['grid_file_path']}", f"{current_working_directory}/grid")
    
    # If it exists, I save it in setup_log and then I copy the grid file from control.dat as "grid"
    #else: 

    #    shutil.copy(f"{current_working_directory}/{p['grid_file_path']}", f"setup_log/grid_old_{dt_string}")
    #    shutil.copy(f"{current_working_directory}/{p['grid_file_path']}", f"{current_working_directory}/grid")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # SOLU FILE
    # if solu file doesn't exist, I copy the solu file from control.dat as "solu"
    #if not os.path.exists(f"{current_working_directory}/solu"):

    shutil.copy(f"{current_working_directory}/{p['solu_file_path']}", f"{current_working_directory}/solu")
    
    # If it exists, I save it in setup_log and then I copy the solu file from control.dat as "solu"
    #else: 

        #shutil.copy(f"{current_working_directory}/{p['solu_file_path']}", f"setup_log/solu_old_{dt_string}")
        #shutil.copy(f"{current_working_directory}/{p['solu_file_path']}", f"{current_working_directory}/solu")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # PHI_INI FILE
    # if phi_ini file doesn't exist, I copy the phi_ini file from control.dat as "phi_ini"
    #if not os.path.exists(f"{current_working_directory}/phi_ini"):

    shutil.copy(f"{current_working_directory}/{p['phi_ini_file_path']}", f"{current_working_directory}/phi_ini")
    
    # If it exists, I save it in setup_log and then I copy the phi_ini file from control.dat as "phi_ini"
    #else: 

    #    shutil.copy(f"{current_working_directory}/{p['phi_ini_file_path']}", f"setup_log/phi_ini_old_{dt_string}")
    #    shutil.copy(f"{current_working_directory}/{p['phi_ini_file_path']}", f"{current_working_directory}/phi_ini")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Check the server and choose the command to execute
    if server_name == 'leftraru2':
        # Server is leftraru2, use sbatch
        execute(f"sbatch {p['queue_filename']}")
    elif server_name == 'Dellemc':
        # Server is Dellemc, use mpirun
        execute(f"mpirun -np {p['number_of_processes']} channel > output.dat 2>&1 &")
    else:
        # Default case if the server name is not recognized
        print(f"Unknown server: {server_name}. Please configure the correct execution command.")

def clean_output(parameters):

    p = parameters
    execute(f"rm -f conver*")
    execute(f"rm -f history*")
    execute(f"rm -f cpu_performance*")
    execute(f"rm -f error.err*")
    execute(f"rm -f run.out*")
    execute(f"rm -f stop.now")
    execute(f"rm -f GlobalMass*")
    execute(f"rm -f GlobalVolume*")
    execute(f"rm -f triangulation*")
    execute(f"rm -f cc*")

    current_directory = os.getcwd()
    directory = f"{current_directory}/{p['output_directory']}"
    os.chdir(directory)

    execute(f"rm -f *.dat")

    execute(f" rm -f combined/*")
    execute(f" rm -f debug/*")
    execute(f" rm -f debugfiles/*")
    execute(f" rm -f solufiles/*")
    execute(f" rm -f solus/*")
    execute(f" rm -f surface/*")
    execute(f" rm -f usolufiles/*")

    os.chdir(current_working_directory)


def clean_setup_log():
    setup_log_directory = f"{current_working_directory}/setup_log"
    os.chdir(setup_log_directory)
    execute(f"rm -f *")
    os.chdir(current_working_directory)


def run_action(action):

    if action == "make_clean":
        make_clean(parameters)
    if action == "compile":
        compile(parameters)
    if action == "setup":

        # If setup_log dir doesn't exist, I create it

        print(f"\n=====================================")
        print(f"Generation of setup files")
        print(f"=====================================\n")

        if not os.path.exists("setup_log"):
            os.makedirs("setup_log")
        control_file_path = actions[1]
        shutil.copy(control_file_path, f"setup_log/{control_file_path.split('.dat')[0]}_{dt_string}.dat")
        write_ind3dmg(parameters)
        write_ind3dmg_LSM(parameters)
        write_directory(parameters)
        write_nproc(parameters)
        write_queue(parameters)
        print(f"\n=====================================")
        print(f"Generation of setup files completed!")
        print(f"=====================================\n")

    if action == "run":
        run(parameters)
    if action == "clean_output":
        clean_output(parameters)
    if action == "clean_setup_log":
        clean_setup_log()


# ------------------
# Main
# ------------------

parameters = get_data()

if "clean_output" in actions:
    run_action("clean_output")
if "clean_setup_log" in actions:
    run_action("clean_setup_log")
if "make_clean" in actions:
    run_action("make_clean")
if "compile" in actions:
    run_action("compile")
if "setup" in actions:
    run_action("setup")
if "run" in actions:
    run_action("run")

