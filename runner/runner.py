import sys, os, json, datetime, os, argparse, resource, subprocess, ntpath, select

RUNNER_DIR = os.path.abspath(os.path.dirname(__file__)) # Directory where runner files are located.
CURRENT_DIR = os.path.abspath(os.getcwd()) # Current directory in the command line.

# Util functions
# Returns: JSON content of the file at the specified path.
def read_json_from_file(file_path):
	with open(file_path, "r") as f:
		return json.loads(f.read())

# Saves the json_object in a file with the specified path.
def save_json_to_file(file_path, json_object):
	with open(file_path, "w") as f:
		f.write(json.dumps(json_object))

# Creates the directory at the specified path if it does not exist.
def create_dir(dir_path):
	if not os.path.isdir(dir_path): os.mkdir(dir_path)

# Returns: text in color blue for console.
def blue(text): return "\033[94m" + text + "\033[0m"
# Returns: text in color red for console.
def red(text): return "\033[91m" + text + "\033[0m"
# Returns: text in color green for console.
def green(text): return "\033[92m" + text + "\033[0m"
# Returns: text in color purple for console.
def purple(text): return "\033[95m" + text + "\033[0m"

def run_program(bin, input_string="", memlim_gb=1024, silent=False):
	# Function that sets limit on memory when executable starts.
	def set_memory_limit():
		soft, hard = resource.getrlimit(resource.RLIMIT_AS)
		resource.setrlimit(resource.RLIMIT_AS, (memlim_gb*1024*1024*1024, hard))

	t_start = datetime.datetime.now()
	process = subprocess.Popen(bin, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE, preexec_fn=set_memory_limit, universal_newlines = True)
	if input_string != "": process.stdin.write(input_string) # Write to STDIN.
	process.stdin.flush()

	stderr_string = ""
	stdout_string = ""
	readable = { process.stdout.fileno(): 1, process.stderr.fileno(): 2 }
	while readable:
		for fd in select.select(readable, [], [])[0]:
			data = os.read(fd, 1024) # read available
			if not data: # EOF
				del readable[fd]
			else:
				if readable[fd] == 1: stdout_string += data.decode("utf-8")
				else:
					stderr_string += data.decode("utf-8")
					if not silent: sys.stdout.buffer.write(data)
					if not silent: sys.stdout.buffer.flush()				
	exit_code = process.wait()
	t_end = datetime.datetime.now()
	process.stdout.close()
	process.stderr.close()

	return {"exit_code": exit_code, "stdout": stdout_string, "stderr": stderr_string, "time": (t_end-t_start).total_seconds() }

# Load runner config file where configurations are stored.
config = read_json_from_file(F"{RUNNER_DIR}/config.json")

# Set directories.
OUTPUT_DIR =  os.path.abspath(os.path.join(RUNNER_DIR,config["output_dir"])) # Directory where output files should be saved.
CMAKELISTS_DIR = os.path.abspath(os.path.join(RUNNER_DIR,config["cmakelists_dir"])) # Directory that contains the root CMakeLists.txt file to compile the project.
INSTANCES_DIR = os.path.abspath(os.path.join(RUNNER_DIR,config["instances_dir"])) # Directory where datasets are stored.
OBJ_DIR = F"{RUNNER_DIR}/obj" # Directory where the object files will be created.

# Set command line parameters.
arg_parser = argparse.ArgumentParser(description="Runs the experiment file(s) specified.")
arg_parser.add_argument("experiments", metavar="EXP_FILE", help="JSON experiment file(s) with the experiments to run.", type=argparse.FileType('r'), nargs='+')
arg_parser.add_argument("--instances", "-I", nargs="*", help="Only execute experiment(s) on selected instances (with these names).")
arg_parser.add_argument("--exps", "-E", nargs="*", help="Only execute selected experiment(s) (with these names).")
arg_parser.add_argument("--callgrind", "-C", help="Runs the experiment(s) using callgrind.", action="store_true")
arg_parser.add_argument("--valgrind", "-V", help="Runs the experiment(s) using valgrind.", action="store_true")
arg_parser.add_argument("--memlimit", "-M", help="Sets a memory limit in GB (default 15GB).", default=15.0)
arg_parser.add_argument("--silent", "-S", help="Do not print the stderr stream of the experiments to the screen.", action="store_true")

# Read command line parameters.
args = vars(arg_parser.parse_args())
experiment_files = args["experiments"]
selected_instances = args["instances"]
selected_experiments = args["exps"]
use_callgrind = args["callgrind"]
use_valgrind = args["valgrind"]
memlimit_gb = float(args["memlimit"])
silent = args["silent"]

# The build type when running callgrind or valgrind is 'debug' otherwise it is 'release'.
build_type = "debug" if use_callgrind or use_valgrind else "release"

# Show parameters.
print(blue(F"Experiment files: {[e.name for e in experiment_files]}"))
print(blue(F"Memory limit: {memlimit_gb}GB"))
print(blue(F"CMakeLists.txt path: {os.path.abspath(CMAKELISTS_DIR)}/CMakeLists.txt"))
print(blue(F"Output directory: {os.path.abspath(OUTPUT_DIR)}"))
print(blue(F"Instances directory: {os.path.abspath(INSTANCES_DIR)}"))

# Returns: the instances that are specified in the experiment file.
# They can be from many datasets and the datasets may have SELECT filters.
# Example: {"name":"ejor2019", "select":"!TAG1|TAG2 TAG3"}
# Meaning: select all instances from dataset ejor2019 that either not containts TAG1 or contains both TAG2 and TAG3.
def instances_for_experiment_file(experiment_file_json):
	instances = []
	for dataset in experiment_file_json["datasets"]:
		dataset_dir = F"{INSTANCES_DIR}/{dataset['name']}"
		dataset_index = read_json_from_file(F"{dataset_dir}/index.json")
		tags_selection_function = lambda tags: True
		if "select" in dataset:
			tag_sets = dataset["select"].split("|")
			tag_sets = [s.split(" ") for s in tag_sets]
			tag_sets = [[t.strip() for t in s if t.strip()] for s in tag_sets if any(t.strip() for t in s)]
			tags_selection_function = lambda tags: any([all([(tag.replace("!", "") in tags) is not ("!" in tag) for tag in tag_set]) for tag_set in tag_sets])
		for entry in dataset_index:
			if selected_instances != None and entry["instance_name"] not in selected_instances: continue # Only leave instances with names matching the filter (if filter is * then all).
			if not tags_selection_function(entry["tags"]): continue # Filter out instances that do not match selection criteria.
			instance = read_json_from_file(F"{dataset_dir}/{entry['file_name']}")
			instance["instance_name"] = entry["instance_name"]
			instance["dataset_name"] = dataset["name"]
			instances.append(instance)
	return instances

# Compiles the CMakeLists.txt in CMAKELISTS_DIR specified in the config.
# Saves the compilation files in OBJ_DIR.
# Returns: if the compilation process was successful.
def compile():
	# Create /obj directory.
	create_dir(OBJ_DIR)
	create_dir(F"{OBJ_DIR}/{build_type}")

	# Compile project using cmake.
	print(purple("Compiling code"), flush=True)
	t0 = datetime.datetime.now()
	os.chdir(F"{OBJ_DIR}/{build_type}")
	exit_code = subprocess.call(["cmake", F"{CMAKELISTS_DIR}", F"-DCMAKE_BUILD_TYPE={build_type}", F"-DRUNNER=ON"])
	if exit_code == 0: exit_code = subprocess.call(["make"])
	os.chdir(CURRENT_DIR)
	if exit_code == 0:
		print(green(F"Finished compiling - Time: {(datetime.datetime.now() - t0).total_seconds()} sec."), flush=True)
	else:
		print(red(F"Compilation failed - Time: {(datetime.datetime.now() - t0).total_seconds()} sec."), flush=True)
	return exit_code == 0

# Run the executable file from the experiment and passes it three arguments to STDIN in the following order:
# experiment: experiment to run.
# instance: instance JSON to use.
# solutions: instance known solutions of the solutions.json file of the dataset.
# Returns: a JSON object with the output of the execution.
# The base format of the JSON object is {"dataset_name":string, "instance_name":string, "experiment_name":string, "success":bool, "stderr":string}.
# 	If success == true: it also has the attributes {"execution_log":json, "solution":json}
#	If success == false: it also has the attributes {"exit_code":number, ""}
def run_experiment(experiment, instance, solutions):
	# Get executable command depending on whether valgrind or callgrind are enabled.
	executable_path = F"{OBJ_DIR}/{build_type}/{experiment['executable']}"
	executable = [executable_path]
	if use_valgrind: executable = ["valgrind", "--track-origins=yes", executable_path]
	elif use_callgrind: executable = ["valgrind", "--tool=callgrind", executable_path]

	# Execute experiment.
	result = run_program(executable, F"{json.dumps(experiment)}{json.dumps(instance)}{json.dumps(solutions)}", memlimit_gb, silent)

	# Try to parse STDOUT as JSON, otherwise leave it as string.
	stdout_json = ""
	try: 
		stdout_json = json.loads(result["stdout"])
	except:
		stdout_json = result["stdout"]

	# If experiment finished successfuly, return the observation and the metadata.
	return {
		"dataset_name":instance["dataset_name"], 
		"instance_name":instance["instance_name"], 
		"experiment_name":experiment["name"], 
		"stderr": result["stderr"], 
		"stdout": stdout_json, 
		"exit_code": result["exit_code"],
		"time": result["time"]
	}

def main():
	# Compile project.
	if not compile(): exit(0)

	# Run experiment files.
	for experiment_file in experiment_files:
		experiment_file_json = json.loads(experiment_file.read())

		# Outputs of the experiments will be stored in this object.
		experiment_file_name = ntpath.basename(experiment_file.name.replace(".json", ""))
		output = {"date": str(datetime.date.today()), "experiment_file": experiment_file_name, "outputs": []}

		# Periodically, every TSave seconds the output will be saved to the output folder with the name "<date>-<experiment_file_name>.json".
		TSave = 60
		output_file_name = F"{datetime.date.today()}-{experiment_file_name}.json"
		TInit = datetime.datetime.now() # TInit = "timestamp when the experimentation started".
		TLast = datetime.datetime.now() # TLast = "last time the output was saved".

		# For each instances specified in the experiment file.
		for instance in instances_for_experiment_file(experiment_file_json):
			# Get instance solutions from the dataset directory.
			solutions = []
			if os.path.isfile(F"{INSTANCES_DIR}/{instance['dataset_name']}/solutions.json"):
				solutions = read_json_from_file(F"{INSTANCES_DIR}/{instance['dataset_name']}/solutions.json")
				solutions = [s for s in solutions if s["instance_name"] == instance["instance_name"]]

			# For each experiment defined in the experiment file.
			for experiment in experiment_file_json["experiments"]:
				# Check if the experiment was selected in the exps argument.
				if selected_experiments != None and experiment["name"] not in selected_experiments: continue

				# Run the experiment.
				print(purple(F"[{instance['dataset_name']}] {instance['instance_name']} - {experiment['name']} ({datetime.datetime.now()})"), flush=True)
				output["outputs"].append(run_experiment(experiment, instance, solutions))
				
				# If TSave seconds have passed since TLast then save output.
				if (datetime.datetime.now() - TLast).total_seconds() >= TSave:
					output["time"] = (datetime.datetime.now() - TInit).total_seconds()
					save_json_to_file(F"{OUTPUT_DIR}/{output_file_name}", output)
					TLast = datetime.datetime.now()
		
		# Having finished all experiments from the experimentation_file, save the final output.
		output["time"] = (datetime.datetime.now() - TInit).total_seconds()
		save_json_to_file(F"{OUTPUT_DIR}/{output_file_name}", output)

if __name__== "__main__":
  main()