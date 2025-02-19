// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// wp_module.cpp: Rcpp R/C++ interface class library -- Rcpp Module example
//
// Copyright (C) 2010 - 2012  Dirk Eddelbuettel and Romain Francois
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>

/// __________________________________________________
/// Count false values in vector<bool> 
int count_falsey(const std::vector<bool> bv)
{
	std::cout << "@count_falsey(const std::vector<bool>)\n";
	auto count = std::count(bv.begin(), bv.end(), false);
    std::cout << "Count = " << count << std::endl;
    return count;
}


/// __________________________________________________
/// Clone vector<double>
std::vector<double> cloney(const std::vector<double>& dv)
{
//	std::cout << "@cloney(const std::vector<double>&) length " << dv.size() << endl;
	std::cout << "@cloney(const std::vector<double>&)\n";
	return std::vector<double>(dv);
}


RCPP_MODULE(wpmod){
    using namespace Rcpp;

    function("count_falsey", &count_falsey, "Count false values in logical vector.");
    function("cloney", &cloney, "Clone vector<double>.");
}
