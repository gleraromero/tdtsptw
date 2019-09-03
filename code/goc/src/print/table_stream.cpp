//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/print/table_stream.h"

#include <math.h>
#include <iomanip>

using namespace std;

namespace goc
{

TableStream::TableStream(ostream* output_stream, double factor)
	:output_stream_(output_stream), factor_(factor), current_attempt_(0), next_valid_attempt_(0)
{ }

TableStream& TableStream::AddColumn(const string& name, int size)
{
	col_names_.push_back(name);
	col_sizes_.push_back(size);
	return *this;
}

bool TableStream::RegisterAttempt()
{
	// If there is no output stream, do nothing.
	return output_stream_ && current_attempt_++ >= next_valid_attempt_;
}

void TableStream::WriteHeader()
{
	if (!output_stream_) return;
	for (int i = 0; i < col_names_.size(); ++i) *output_stream_ << setw(col_sizes_[i]) << col_names_[i];
	*output_stream_ << endl;
}

void TableStream::WriteRow(const vector<string>& row)
{
	// If there is no output stream, do nothing.
	if (!output_stream_) return;
	
	// Output current row.
	for (int i = 0; i < row.size(); ++i) *output_stream_ << setw(col_sizes_[i]) << row[i];
	if (output_stream_) *output_stream_ << endl;
	
	// Increment next valid attempt by the factor.
	next_valid_attempt_ = (int)floor((double)current_attempt_ * factor_);
}
} // namespace goc