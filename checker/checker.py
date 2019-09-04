import sys, os, json, datetime, os, argparse

CHECKER_DIR = os.path.abspath(os.path.dirname(__file__)) # Directory where checker files are located.
CURRENT_DIR = os.path.abspath(os.getcwd()) # Current directory in the command line.
INSTANCES_DIR = os.path.abspath(os.path.join(CHECKER_DIR, "..", "instances")) # Directory where datasets are stored.

# Util functions
# Returns: JSON content of the file at the specified path.
def read_json_from_file(file_path):
	with open(file_path, "r") as f:
		return json.loads(f.read())

# Returns: text in color blue for console.
def blue(text): return F"\033[94m{text}\033[0m"
# Returns: text in color red for console.
def red(text): return F"\033[91m{text}\033[0m"
# Returns: text in color green for console.
def green(text): return F"\033[92m{text}\033[0m"
# Returns: text in color purple for console.
def purple(text): return F"\033[95m{text}\033[0m"

INFTY = 10e8; EPS = 10e-4
def epsilon_equal(a, b): return abs(a - b) < EPS
def epsilon_bigger(a, b): return a > b + EPS

# Set command line parameters.
arg_parser = argparse.ArgumentParser(description="Check if the solutions in the output file are correct.")
arg_parser.add_argument("output_files", metavar="OUT_FILE", help="JSON output file(s) to run the checker on.", type=argparse.FileType('r'), nargs='+')

# Read command line parameters.
args = vars(arg_parser.parse_args())
output_files = args["output_files"]

# Show parameters.
print(blue(F"Output files: {[e.name for e in output_files]}"))

file_name_by_instance = {}

def read_instance(dataset_name, instance_name):
	if not dataset_name in file_name_by_instance:
		index = read_json_from_file(F"{INSTANCES_DIR}/{dataset_name}/index.json")
		file_name_by_instance[dataset_name] = {}
		for entry in index:
			file_name_by_instance[dataset_name][entry["instance_name"]] = F"{INSTANCES_DIR}/{dataset_name}/{entry['file_name']}"
	return read_json_from_file(file_name_by_instance[dataset_name][instance_name])

# Returns: the best known solution (with minimum value) among all with the specific tags.
solutions_cache = {}
def best_known_solution(dataset_name, instance_name, tags):
	if dataset_name not in solutions_cache:
		solutions_cache[dataset_name] = read_json_from_file(F"{INSTANCES_DIR}/{dataset_name}/solutions.json")
	solutions = solutions_cache[dataset_name]
	bks = {"value":10e8}
	for s in solutions:
		if s["instance_name"] == instance_name and s["tags"] == tags and s["value"] <= bks["value"]:
			bks = s
	return bks

# Returns: the travel time of arc (i, j) if departing at t0.
def travel_time(instance, i, j, t0):
    c = instance["clusters"][i][j]
    T = instance["speed_zones"]
    v = instance["cluster_speeds"]

    # Find speed slot T_k that includes t0.
    k = min(i for i in range(0,len(T)) if t0 >= T[i][0] and t0 <= T[i][1])

    # Travel time algorithm from Ichoua et al.
    t = t0
    d = instance["distances"][i][j]
    tt = t + d / v[c][k]
    while tt > T[k][1]:
        d = d - v[c][k] * (T[k][1] - t)
        t = T[k][1]
        if epsilon_equal(d, 0): break
        if epsilon_bigger(d, 0) and k+1 == len(T): return INFTY
        tt = t + d / v[c][k+1]
        k = k + 1
    return tt - t0

# Returns: the time when the path finishes if departing at t0.
# If the route is infeasible it returns INFTY.
def ready_time(instance, path, t0):
	n = instance["digraph"]["vertex_count"] # Number of vertices.
	T = instance["horizon"][1] # End of horizon.
	tw = instance["time_windows"] # Time windows
	a = [tw[i][0] for i in range(0, n)] # Release times
	b = [tw[i][1] for i in range(0, n)] # Deadline times

	# Check path elementarity.
	if len(path) != len(set(path)): return INFTY

	# Check that all vertices are visited if instance is not profitable.
	if len(path) != n: return INFTY

	# Calculate ready time.
	t = max(t0, a[path[0]])
	for k in range(0, len(path)-1):
		i = path[k] 
		j = path[k+1]
		t = max(t+travel_time(instance, i, j, t), a[j])
		if epsilon_bigger(t, b[j]): return INFTY
	return t

# Returns: true if the routes are valid.
# If routes are not valid, then it adds the error messages to error_messages
def check_routes(instance, solution, error_messages):
	valid = True
	for route in solution["routes"]:
		tf = ready_time(instance, route["path"], route["t0"])
		expected_duration = tf - route["t0"]
		if tf == INFTY: 
			error_messages.append(F"\tInfeasible: {route}")
			valid = False
		elif not epsilon_equal(expected_duration, route["duration"]): 
			error_messages.append(F"\tDifferent duration: {route} - Expected: {expected_duration}")
			valid = False
	return valid

def main():
	ok = 0
	wrong = 0
	suboptimal = 0
	skip = 0
	mlim = 0
	err = 0
	for output_file in output_files:
		output_file_json = json.loads(output_file.read())
		print(F"Checking {output_file.name}")

		for output in output_file_json["outputs"]:
			dataset_name = output["dataset_name"]
			instance_name = output["instance_name"]
			experiment_name = output["experiment_name"]
			if output["exit_code"] == -6:
				mlim += 1
				print(F"Checking {experiment_name} - {dataset_name} {instance_name}: {purple('Memory limit.')}")
				continue
			if output["exit_code"] != 0:
				err += 1
				print(F"Checking {experiment_name} - {dataset_name} {instance_name}: {purple('Exit code: ' + str(output['exit_code']))}")
				continue

			instance = read_instance(dataset_name, instance_name)
			tags = []
			if "makespan" in experiment_name: tags.append("MAKESPAN")
			if "duration" in experiment_name: tags.append("DURATION")
			bks = best_known_solution(dataset_name, instance_name, tags)

			status = output["stdout"]["Exact"]["status"]
			opt_found = status == "Optimum" or status == "Finished" # Indicates if the optimum solution was found.
			if "Best solution" in output["stdout"]:
				solution = output["stdout"]["Best solution"]
				errors = []
				valid = check_routes(instance, solution, errors)
				if not valid: 
					print(F"Checking {experiment_name} - {dataset_name} {instance_name}")
					for error in errors: print(red(error))
					wrong += 1
				elif opt_found and epsilon_bigger(solution["value"], bks["value"]):
					suboptimal += 1
					print(F"Checking {experiment_name} - {dataset_name} {instance_name}")
					print(blue(F"Suboptimal - BKS: {bks} - Obtained: {solution}"))
				else:
					ok += 1
			else:
				print(F"Checking {experiment_name} - {dataset_name} {instance_name}: {purple('No solution (' + status + ')')}")
				skip += 1

	print(green(F"ok: {ok}"), red(F"wrong: {wrong}"), blue(F"suboptimal: {suboptimal}"), purple(F"skipped: {skip}"), purple(F"memlim: {mlim}"), purple(F"error: {err}"))

if __name__== "__main__":
  main()