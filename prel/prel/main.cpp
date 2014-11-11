//
//  main.cpp
//  prel
//
//  Created by George Oblapenko on 11/11/14.
//  Copyright (c) 2014 George Oblapenko. All rights reserved.
//

#include <iostream>
#include <kineticlib.h>

int main(int argc, const char * argv[]) {
    double start_T = 2000.;  // the start temperature
    double end_T = 10000.0;  // the end temperature
    int points_amt = 50;  // the amount of points
    double xN = 0.2;  // the relative numeric density of atomic nitrogen
    
    arma::vec T_array = arma::linspace<arma::vec>(start_T, end_T, points_amt);  // start value of T, end value of T, amount of steps
    arma::vec prel = arma::zeros(50);  // array of relaxation pressure values
    double p = 100000.0;  // atmospheric pressure
    
    klib::MoleculeOneT N2 = klib::MoleculeOneT("N2");
    klib::Atom N = klib::Atom("N");
    
    arma::vec idata_N2N2 = klib::load_elastic_parameters(N2, N2);  // N2 + N2 elastic interaction data
    arma::vec idata_N2N = klib::load_elastic_parameters(N2, N);  // N2 + N elastic interaction data
    
    arma::vec ddata_N2N2 = klib::load_dissociation_parameters(N2, N2);  // N2 + N2 dissociation data
    arma::vec ddata_N2N = klib::load_dissociation_parameters(N2, N);  // N2 + N dissociation data
    
    arma::mat33 beta_matrix;  // the matrix of beta integral brackets
    arma::vec3 right_parts;  // the right-hand side of the system
    arma::vec3 results;  // the values of the expansion coefficients g_{c,pq} - the order is g_{N2,10}, g_{N2,01}, g_{N,10}
    
    double T;
    double n;  // the numeric density of the mixture
    double rho;  // the density of the mixture
    
    beta_matrix.at(0, 0) = 1.5 * (1. - xN);
    beta_matrix.at(0, 2) = 1.5 * xN;  // these appear due to the constraint conditions and are independent of temperature
    
    for (int i=0; i<points_amt; i++) {
        T = T_array[i];
        n = p / T;
        N2.renorm(T, (1. -xN) * n);
        rho = (1. -xN) * n * N2.mass + xN * n * N.mass;
        
    }
}
