//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include <iostream>
#include <vector>
#include <goc/goc.h>

#include "vrp_instance.h"
#include "preprocess_travel_times.h"
#include "preprocess_time_windows.h"
#include "preprocess_waiting_times.h"
#include "labeling.h"

using namespace std;
using namespace goc;
using namespace nlohmann;
using namespace tdtsptw;

int main(int argc, char** argv)
{
	json output; // STDOUT output will go into this JSON.
	
	simulate_runner_input("instances/guerriero_et_al_2014b", "40_70_A_100_A1", "experiments/easy.json", "makespan");
	
	json experiment, instance, solutions;
	cin >> experiment >> instance >> solutions;
	
	// Parse experiment.
	Duration time_limit = Duration(value_or_default(experiment, "time_limit", 7200), DurationUnit::Seconds);
	string objective = value_or_default(experiment, "objective", "makespan");
	
	// Show experiment details.
	clog << "Time limit: " << time_limit << "sec." << endl;
	clog << "Objective: " << objective << endl;
	
	// Preprocess instance JSON.
	preprocess_travel_times(instance);
	preprocess_waiting_times(instance);
	preprocess_time_windows(instance);
	
	// Parse instance.
	VRPInstance vrp = instance;
	MLBExecutionLog log;
	Route best = run_labeling(vrp, time_limit, &log, objective == "makespan");
	
	// Show results.
	clog << "Time: " << log.time << endl;
	clog << "Status: " << log.status << endl;
	if (!best.path.empty())
	{
		clog << "Best route: " << best.path << ", departing at: " << best.t0 << ", duration: " << best.duration << endl;
		output["Best solution"] = VRPSolution(best.duration, {best});
	}
	output["Exact"] = log;

	// Send JSON output to cout.
	cout << output << endl;
	return 0;
}