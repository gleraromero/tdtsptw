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

using namespace std;
using namespace goc;
using namespace nlohmann;
using namespace tdtsptw;

int main(int argc, char** argv)
{
	json output; // STDOUT output will go into this JSON.
	
	simulate_runner_input("instances/guerriero_et_al_2014", "15_70_A_A1", "experiments/tdtsptw.json", "Main");
	
	json experiment, instance, solutions;
	cin >> experiment >> instance >> solutions;
	
	// Parse experiment.
	Duration time_limit = Duration(value_or_default(experiment, "time_limit", 7200), DurationUnit::Seconds);
	
	// Show experiment details.
	clog << "Time limit: " << time_limit << "sec." << endl;
	
	// Preprocess instance JSON.
	preprocess_travel_times(instance);
	preprocess_time_windows(instance);
	
	// Parse instance.
	VRPInstance vrp = instance;
	MLBExecutionLog log(true);
	Route best;
	
	// Show results.
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