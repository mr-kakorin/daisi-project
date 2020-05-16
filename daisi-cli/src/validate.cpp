//
// Created by artoria on 3/15/20.
//

#include "validate.h"

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>

void validate( const std::vector<int> &middles, const std::vector<int> &lefts, const std::vector<int> &rights,
               const std::vector<int> &downs, const std::vector<int> &ups,
               const std::vector<double> &c_middles, const std::vector<double> &c_lefts,
               const std::vector<double> &c_rights, const std::vector<double> &c_downs,
               const std::vector<double> &c_ups,
               const std::vector<double> &vect,
               std::vector<double> &x,
               const std::vector<double> &x_valid ) {
	short flag = 0;
	double diff_sum = 0, tmp_sum = 0;
	double w = 1.9;
	double eps_tolerance = 1e-9;
	int n = 0;

	while (flag == 0) {
		diff_sum = 0;

		for (int k = 0; k < middles.size(); ++k) {
			n = middles[k];
			tmp_sum = -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] +
			            x[ups[k]] * c_ups[k] + x[downs[k]] * c_downs[k]) *
			          w / c_middles[k];
			tmp_sum = tmp_sum + (1 - w) * x[n] + w * vect[n] / c_middles[k];
			if ((x[n] > eps_tolerance * 0.1) || (x[n] < -eps_tolerance * 0.1)) {
				diff_sum = diff_sum + std::abs( x[n] - tmp_sum ) / x[n];
			} else {
				diff_sum = diff_sum + std::abs( x[n] - tmp_sum );
			}
			x[n] = tmp_sum;
		}
		//progress = std::abs( diff_sum ) / eps_tolerance;
		if (std::abs( diff_sum ) < eps_tolerance)
			flag = 1;
	}

	for (int i = 0; i < x.size(); ++i) {
		assert( std::abs( x[i] - x_valid[i] ) < 1e-9 );
	}
}

void validate_parallel( const std::vector<int> &middles, const std::vector<int> &lefts, const std::vector<int> &rights,
                        const std::vector<int> &downs, const std::vector<int> &ups,
                        const std::vector<double> &c_middles, const std::vector<double> &c_lefts,
                        const std::vector<double> &c_rights, const std::vector<double> &c_downs,
                        const std::vector<double> &c_ups,
                        const std::vector<double> &vect,
                        std::vector<double> &x,
                        const std::vector<double> &x_valid ) {
	short flag = 0;
	double diff_sum = 0, tmp_sum = 0;
	double w = 1.9;
	double eps_tolerance = 1e-9;
	int n = 0;
	std::vector<double> z;
	z.resize( x.size());
	while (flag == 0) {
		diff_sum = 0;

		for (int k = 0; k < middles.size(); ++k) {
			n = middles[k];
			z[n] = (1 - w) * x[n];
			z[n] -= ( w / c_middles[k] * c_ups[k]) * x[ups[k]];
			z[n] -= ( w / c_middles[k] * c_rights[k]) * x[rights[k]];
			z[n] += w * vect[n] / c_middles[k];
		}

		for (int k = 0; k < middles.size(); ++k) {
			n = middles[k];
			tmp_sum = z[n] - (1 + w / c_middles[k] * (x[lefts[k]] * c_lefts[k] + x[downs[k]] * c_downs[k]));

			if ((x[n] > eps_tolerance * 0.1) || (x[n] < -eps_tolerance * 0.1)) {
				diff_sum = diff_sum + std::abs( x[n] - tmp_sum ) / x[n];
			} else {
				diff_sum = diff_sum + std::abs( x[n] - tmp_sum );
			}
			x[n] = tmp_sum;
		}


//		for (int k = 0; k < middles.size(); ++k)
//		{
//			n = middles[k];
//			tmp_sum = -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] +
//			            x[ups[k]] * c_ups[k] + x[downs[k]] * c_downs[k]) *
//			          w / c_middles[k];
//			tmp_sum = tmp_sum + (1 - w) * x[n] + w * vect[n] / c_middles[k];
//			if ((x[n] > eps_tolerance * 0.1) || (x[n] < -eps_tolerance * 0.1)) {
//				diff_sum = diff_sum + std::abs( x[n] - tmp_sum ) / x[n];
//			} else {
//				diff_sum = diff_sum + std::abs( x[n] - tmp_sum );
//			}
//			x[n] = tmp_sum;
//		}
		//progress = std::abs( diff_sum ) / eps_tolerance;
		if (std::abs( diff_sum ) < eps_tolerance)
			flag = 1;
	}

	for (int i = 0; i < x.size(); ++i) {
		assert( std::abs( x[i] - x_valid[i] ) < 1e-9 );
	}
}


template<typename T>
void read_vec_from_file(const std::string& in_path, std::vector<T>& out_vec)
{
	std::ifstream ifile(in_path, std::ios::in);
	//check to see that the file was opened correctly:
	if (!ifile.is_open()) {
		std::cerr << "There was a problem opening the input file!\n";
		exit(1);//exit or do additional error checking
	}
	T num;
	while (ifile >> num) {
		out_vec.push_back(num);
	}
}


void run()
{
	std::vector<int> middles, lefts, rights, downs, ups;
	std::vector<double> c_middles, c_lefts, c_rights, c_downs, c_ups, vect, x, x_valid;
	std::string base_path( "/home/artoria/projects/daisi/daisi-cli/validate/" );

	read_vec_from_file<int>( base_path + "middles.csv", middles );
	read_vec_from_file<double>( base_path + "c_middles.csv", c_middles );

	read_vec_from_file<int>( base_path + "lefts.csv", lefts );
	read_vec_from_file<double>( base_path + "c_lefts.csv", c_lefts );

	read_vec_from_file<int>( base_path + "rights.csv", rights );
	read_vec_from_file<double>( base_path + "c_rights.csv", c_rights );

	read_vec_from_file<int>( base_path + "downs.csv", downs );
	read_vec_from_file<double>( base_path + "c_downs.csv", c_downs );

	read_vec_from_file<int>( base_path + "ups.csv", ups );
	read_vec_from_file<double>( base_path + "c_ups.csv", c_ups );

	read_vec_from_file<double>( base_path + "vect.csv", vect );
	read_vec_from_file<double>( base_path + "x.csv", x_valid );
	read_vec_from_file<double>( base_path + "x_begin.csv", x );

	validate_parallel( middles, lefts, rights, downs, ups,
	                   c_middles, c_lefts, c_rights, c_downs, c_ups,
	                   vect, x, x_valid );
}