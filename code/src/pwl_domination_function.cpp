//
// Created by gleraromero on 27/06/19.
//

#include "pwl_domination_function.h"

using namespace std;
using namespace goc;

namespace tdtsptw
{
PWLDominationFunction::PWLDominationFunction(const PWLFunction& f)
{
	// Add all f's pieces in order.
	size_ = 0;
	int i = first_ = last_ = -1;
	pieces_.reserve(f.Pieces().size()*2);
	next_.reserve(f.Pieces().size()*2);
	for (auto& p: f.Pieces()) i = AddPieceAfter(i, p);
	domain_ = f.Domain();
}

PWLDominationFunction::operator PWLFunction() const
{
	PWLFunction f;
	for (int i = first_; i != -1; i = next_[i]) f.AddPiece(pieces_[i]);
	return f;
}

Interval PWLDominationFunction::Domain() const
{
	return domain_;
}

bool PWLDominationFunction::Empty() const
{
	return size_ == 0;
}

bool PWLDominationFunction::DominatePieces(const PWLFunction& f2, double delta)
{
	if (Empty()) return true;
	auto& f1 = *this;
	if (epsilon_smaller(f1.Domain().right, f2.Domain().left)) return false;
	
	// Check domination if functions are points.
	if (f2.Domain().IsPoint() && !f1.Domain().IsPoint()) return false;

	// Piece to add to the final of f2 with waiting time to include all f1's domain.
	double f2_last_duration = f2.LastPiece().Value(f2.Domain().right);
	auto completion_piece = LinearFunction(
		{max(dom(f2)), f2_last_duration},
		{f1.Domain().right, f2_last_duration + (f1.Domain().right - max(dom(f2)))}
	);
	
	int i1 = first_, i2 = 0;
	int prev_i1 = -1; // index of the previous element of i1.
	while (i1 != -1 && i2 <= f2.PieceCount())
	{
		auto& p1 = f1.pieces_[i1];
		auto& p2 = i2 == f2.PieceCount() ? completion_piece : f2.Piece(i2);
		
		// If p2 = [ .... ] --nointersection- [ ....] = p1, move p2 forward.
		if (epsilon_smaller(max(dom(p2)), min(dom(p1)))) { ++i2; continue; }
		
		// If p1 = [ .... ] --nointersection- [ .... ] = p2, move p1 forward.
		if (epsilon_smaller(p1.domain.right, p2.domain.left)) { prev_i1 = i1; i1 = next_[i1]; continue; }
		
		// If p2 intersects p1 exactly at the end, and p1 is not a point then it makes
		// no sense in attempting domination.
		if (epsilon_equal(p1.domain.right, p2.domain.left) && !p1.domain.IsPoint()) { prev_i1 = i1; i1 = next_[i1]; continue; }
		
		// If p2 intersects p1 exactly at the begining, and p1 is not a point then it makes
		// no sense in attempting domination.
		if (epsilon_equal(p2.domain.left, p1.domain.right) && !p1.domain.IsPoint()) { ++i2; continue; }
		
		// Now we have intersection between p1 and p2, we call it [l, r].
		double l = max(p1.domain.left, p2.domain.left);
		double r = min(p1.domain.right, p2.domain.right);
		
		double f1l = p1.Value(l), f1r = p1.Value(r);
		double f2l = p2.Value(l), f2r = p2.Value(r);
		
		// There is some domination if f1(l) >= f2(l)+delta or f1(r) >= f2(r)+delta.
		if (epsilon_bigger_equal(f1l, f2l+delta) || epsilon_bigger_equal(f1r, f2r+delta))
		{
			// Save the dominated interval in [ld, rd].
			double ld = l, rd = r;
			// Case 1: f1 is dominated in all [l, r].
			if (epsilon_bigger_equal(f1l, f2l+delta) && epsilon_bigger_equal(f1r, f2r+delta))
			{
				ld = l;
				rd = r;
			}
				// Case 2: f1 is dominated in a prefix or suffix of [l, r].
			else
			{
				// Find intersection in domain m such that f1(m) == f2(m)+delta.
				double m = p1.Intersection(LinearFunction({l, f2l+delta}, {r, f2r+delta}));
				m = max(m, l);
				m = min(m, r);
				// If prefix is dominated then f1(l) >= f1(l)+delta
				if (epsilon_bigger_equal(f1l, f2l+delta)) rd = m;
					// If suffix is dominated then f1(r) >= f2(r)+delta.
				else if (epsilon_bigger_equal(f1r, f2r+delta)) ld = m;
			}
			
			// Dominate p1 in interval [ld, rd].
			// Case A: [ld, rd] is a superset of dom(p1), then remove p1.
			if (epsilon_smaller_equal(ld, p1.domain.left) && epsilon_bigger_equal(rd, p1.domain.right))
			{
				i1 = ErasePiece(i1, prev_i1);
				// We moved i1 to the next piece, we can move to the next iteration.
				continue;
			}
				// Case B: [ld, rd] is a prefix of dom(p1), then update the left domain and image of p1.
			else if (epsilon_smaller_equal(ld, p1.domain.left))
			{
				p1 = LinearFunction({rd, p1.Value(rd)}, {p1.domain.right, p1.Value(p1.domain.right)});
			}
				// Case C: [ld, rd] is a suffix of dom(p1), then update the right domain and image of p1.
			else if (epsilon_bigger_equal(rd, p1.domain.right))
			{
				p1 = LinearFunction({p1.domain.left, p1.Value(p1.domain.left)}, {ld, p1.Value(ld)});
			}
				// Case D: [ld, rd] is in the middle of dom(p1), then we need to split p1 into the leftmost and rightmost pieces.
				// Only apply this case if ld < rd. Discarding a point in a continuous space is unnecessary.
			else if (epsilon_smaller(ld, rd))
			{
				// We first create the rightmost piece [rd, max(dom(p1))].
				// Move i1 to this piece.
				prev_i1 = i1;
				auto right_piece = LinearFunction({rd, p1.Value(rd)}, {p1.domain.right, p1.image.right});
				i1 = AddPieceAfter(i1, right_piece);
				
				// Now we shrink the leftmost piece to [min(dom(p1)), ld].
				pieces_[prev_i1] = LinearFunction({pieces_[prev_i1].domain.left, pieces_[prev_i1].Value(pieces_[prev_i1].domain.left)}, {ld, pieces_[prev_i1].Value(ld)});
			}
		}
		
		// Domination is done, now move the piece which ends before.
		if (pieces_[i1].domain.right < p2.domain.right) { prev_i1 = i1; i1 = next_[i1]; }
		else i2++;
	}
	return Empty();
}

bool PWLDominationFunction::IsAlwaysDominated(const PWLFunction& f2, double delta)
{
	if (Empty()) return true;
	PWLDominationFunction& f1 = *this;
	
	// If earliest arrival of f2 is later than f1's then all domain can not be dominated.
	if (epsilon_bigger(min(dom(f2)), f1.Domain().left)) return false;
	
	// Piece to add to the final of f2 with waiting time to include all f1's domain.
	double f2_last_duration = f2.LastPiece().Value(f2.Domain().right);
	auto completion_piece = LinearFunction(
		{max(dom(f2)), f2_last_duration},
		{f1.Domain().right, f2_last_duration + (f1.Domain().right - max(dom(f2)))}
	);
	
	int i1 = first_, i2 = 0;
	while (i1 != -1 && i2 <= f2.Pieces().size())
	{
		auto& p1 = f1.pieces_[i1];
		
		// Search first piece in f2 that ends after the begining of p1.
		for (; i2 < f2.Pieces().size(); ++i2)
			if (epsilon_bigger_equal(f2.Piece(i2).domain.right, p1.domain.right))
				break;
		auto& p2 = i2 == f2.PieceCount() ? completion_piece : f2.Piece(i2);
		
		// If p1 = [ .... ] --space-- [ .... ] = p2, then p1 is not dominated.
		if (epsilon_smaller(p1.domain.right, p2.domain.left)) return false;
		
		// Now we have intersection between p1 and p2, we call it [l, r].
		double l = max(p1.domain.left, p2.domain.left);
		double r = min(p1.domain.right, p2.domain.right);
		
		double f1l = p1.Value(l), f1r = p1.Value(r);
		double f2l = p2.Value(l), f2r = p2.Value(r);
		
		// If f1(l) < f2(l)+delta or f1(r) < f2(r)+delta, then there is some part of f1 which is not dominated.
		if (epsilon_smaller(f1l, f2l+delta) || epsilon_smaller(f1r, f2r+delta)) return false;
		
		// Move the piece which ends before.
		if (p1.domain.right < p2.domain.right) i1 = next_[i1];
		else i2++;
	}
	return true;
}

int PWLDominationFunction::ErasePiece(int i, int prev_i)
{
	if (i == last_) last_ = prev_i;
	if (prev_i != -1) next_[prev_i] = next_[i];
	else first_ = next_[i];
	size_--;
	if (size_ == 0) domain_ = {INFTY, -INFTY};
	else domain_ = {pieces_[first_].domain.left, pieces_[last_].domain.right};
	return next_[i];
}

int PWLDominationFunction::AddPieceAfter(int i, const LinearFunction& piece)
{
	int j = pieces_.size();
	pieces_.push_back(piece);
	next_.push_back(-1);
	if (i != -1)
	{
		next_[j] = next_[i];
		next_[i] = j;
	}
	if (i == last_) last_ = j;
	if (i == -1) first_ = j;
	size_++;
	domain_ = {pieces_[first_].domain.left, pieces_[last_].domain.right};
	return j;
}
} // namespace tdtsptw