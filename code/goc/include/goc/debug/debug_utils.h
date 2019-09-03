//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_DEBUG_DEBUG_UTILS_H
#define GOC_DEBUG_DEBUG_UTILS_H

#include <string>

namespace goc
{
// This is a function that is useful if debugging in an IDE which does not allow to specify an input for the STDIN.
// This function simulates that the program was called by the runner.py script with the instance 'instance_name' in
// the dataset in the directory 'dataset_dir'. Also it receives the experiment file located in the 'experiment_path' and
// inputs the 'experiment_name'. Also if the instance has solutions in the corresponding solutions.json file then it
// inputs those solutions too.
// Effect: adds to STDIN << experiment << instance << solutions;
void simulate_input_in_debug(const std::string& dataset_dir, const std::string& instance_name,
							 const std::string& experiment_path, const std::string& experiment_name);
} // namespace goc

#endif // GOC_DEBUG_DEBUG_UTILS_H