//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_PRINT_TABLE_STREAM_H
#define GOC_PRINT_TABLE_STREAM_H

#include <iostream>
#include <string>
#include <vector>

namespace goc
{
// This is a stream for outputting tables with information. The idea is that not all rows should be displayed, and only
// a subset of them to not overflow the screen. This is for outputting execution information as an algorithm with
// multiple iterations runs.
//
// Example:
// 	AddColumn("Name", 10);
// 	AddColumn("Age", 5);
// 	WriteHeader();
// 	WriteRow("Mike", 21);
// 	WriteRow("Jessica", 43);
// Output:
//	Name      Age
//	Mike	  21
//	Jessica	  43
class TableStream
{
public:
	// output_stream: ostream where the table will be outputted.
	// factor: rows will be allowed to be written if their number is a power of 'factor'.
	TableStream(std::ostream* output_stream, double factor);
	
	// Add column to the table stream.
	// Returns: a reference to this object to concatenate multiple calls.
	TableStream& AddColumn(const std::string& name, int size);
	
	// Registers the attempt to write a row.
	// Returns: if this attempt should be actually written.
	bool RegisterAttempt();
	
	// Writes the header of the table.
	void WriteHeader();
	
	// Writes a row to the table with the values in the row associated to each column in order.
	void WriteRow(const std::vector<std::string>& row);

private:
	std::ostream* output_stream_; // the table is outputted here.
	double factor_; // an attempt to write a row is valid if it is a power of 'factor'.
	int current_attempt_; // the current attempt number to write a row.
	int next_valid_attempt_; // the next attempt that will be valid to be written.
	std::vector<std::string> col_names_; // name of the columns.
	std::vector<int> col_sizes_; // sizes of the columns.
};
} // namespace goc

#endif //GOC_PRINT_TABLE_STREAM_H
