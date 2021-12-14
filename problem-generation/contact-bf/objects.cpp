/*

    This file contains the basic
    objects employed by the rest 
    of the code. 

*/

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <assert.h>
#include <string>

#include "objects.h"


/*
    Sequence methods
*/

Sequence::Sequence(int len, int dim) : len(len), dim(dim) {

    assert(len > 0);

}

void Sequence::set_seq(unsigned long long int k) {
    /*
        This method converts a binary number to
        base 4 (or 6, if it is a 3D sequence), and 
        maps every number to a string of turns.
    */

    assert(seq.empty());
    assert(dim == 2 or dim == 3);

    unsigned long long t = k;

    if (dim == 2) {

        for (int i=0; i<len; i++) {

            int rem = t % 4;
            switch (rem) {

                case 0:
                    seq.push_back('U');
                    break;
                case 1:
                    seq.push_back('D');
                    break;
                case 2:
                    seq.push_back('R');
                    break;
                case 3:
                    seq.push_back('L');
                    break;

                default:
                    std::cout << "Fatal error: mod 4 operation evaluated " << rem << ".\n";
                    exit(1);

            }

            t = (unsigned long long) floor(t / 4);

        }

    } else if (dim == 3) {

        for (int i=0; i<len; i++) {

            int rem = t % 6;
            switch (rem) {

                case 0:
                    seq.push_back('U');
                    break;
                case 1:
                    seq.push_back('D');
                    break;
                case 2:
                    seq.push_back('R');
                    break;
                case 3:
                    seq.push_back('L');
                    break;
                case 4:
                    seq.push_back('I');
                    break;
                case 5:
                    seq.push_back('O');
                    break;

                default:
                    std::cout << "Fatal error: mod 6 operation evaluated " << rem << ".\n";
                    exit(1);

            }

            t = (unsigned long long) floor(t / 6);

        }

    } else {

        std::cout << "Fatal error: dimension is not 2 or 3.\n";
        std::cout << "(and this was uncaught by assertion)\n";
        exit(1);

    }

    std::cout << k << " : " << seq << "\n";
    t = k;
    for (int i=0; i<len; i++) {
        int rem = t % 6;
        std::cout << rem << " ";
        t = (unsigned long long) floor(t / 6);
    }
    std::cout << "\n";

}

int Sequence::get_len() {

    return len;

}

int Sequence::get_dim() {

    return dim;

}

std::string Sequence::get_seq() {

    return seq;

}

void Sequence::reset() {

    seq.clear();

}

/*
    ContactSet methods
*/

ContactSet::ContactSet(std::string id) : id(id) {

    n_contacts = 0;

}

void ContactSet::add_contact(int res1, int res2) {

    first.push_back(res1);
    second.push_back(res2);

    n_contacts++;

}

int ContactSet::get_n_contacts() {

    return n_contacts;

}

std::string ContactSet::get_id() {

    assert(!id.empty());
    return id;

}

std::list <int> ContactSet::get_first_list() {

    return first;

}

std::list <int> ContactSet::get_second_list() {

    return second;

}

/*
    Peptide methods
*/

Peptide::Peptide(int n_res, int dim) : n_res(n_res), dim(dim) {

    assert(dim == 2 or dim == 3);

    if (dim == 2) {

        x = new int[n_res];
        y = new int[n_res];

        x[0] = 0;
        y[0] = 0;
        z = NULL;
        res_placed = 1;

    } else if (dim == 3) {

        x = new int[n_res];
        y = new int[n_res];
        z = new int[n_res];

        x[0] = 0;
        y[0] = 0;
        z[0] = 0;
        res_placed = 1;

    } else {

        std::cout << "Fatal error: dimension is not 2 or 3.\n";
        std::cout << "(and this was uncaught by assertion)\n";
        exit(1);

    }
}

Peptide::~Peptide() {


    if (dim == 2) {

        assert(z == NULL);

        delete[] x;
        delete[] y;

    } else if (dim == 3) {

        delete[] x;
        delete[] y;
        delete[] z;

    } else {

        std::cout << "Fatal error: dimension is not 2 or 3.\n";
        exit(1);

    }

}

int Peptide::last_x() {

    return x[res_placed-1];

}

int Peptide::last_y() {

    return y[res_placed-1];

}

int Peptide::last_z() {

    assert(dim == 3);
    return z[res_placed-1];

}

int Peptide::get_dim() {

    return dim;

}

int Peptide::get_seq_length() {

    return n_res;

}

int Peptide::get_chain_length() {

    return res_placed;

}

bool Peptide::is_clashing(int _x, int _y, int _z) {

    if (dim == 2) {

        for (int i=0; i<res_placed; i++) {
            if (x[i] == _x && y[i] == _y) {
                return true;
            }
        }
        return false;

    } else if (dim == 3) {

        for (int i=0; i<res_placed; i++) {
            if (x[i] == _x && y[i] == _y && z[i] == _z) {
                return true;
            }
        }
        return false;

    } else {

        std::cout << "Fatal error: dimension is not 2 or 3.\n";
        exit(1);

    }

}

void Peptide::add_residue(int _x, int _y, int _z) {

    if (n_res == res_placed) {
        std::cout << "Attempted to add residues to finished SAW.\n";
        exit(1);
    }

    if (dim == 2) {
        x[res_placed] = _x;
        y[res_placed] = _y;
        res_placed++;

    } else if (dim == 3) {

        x[res_placed] = _x;
        y[res_placed] = _y;
        z[res_placed] = _z;
        res_placed++;

    } else {

        std::cout << "Fatal error: dimension is not 2 or 3.\n";
        exit(1);

    }

}

ContactSet Peptide::get_contact_set(std::string seq) {

    ContactSet set(seq);

    if (dim == 2) {

        for (int i=0; i<n_res; i++) {
            for (int j=i+3; j<n_res; j++) {

                if (abs(j-i) % 2 == 0 or (j-1) == 1) continue;

                if ((abs(x[i]-x[j]) == 1 && abs(y[i]-y[j]) == 0) xor
                    (abs(x[i]-x[j]) == 0 && abs(y[i]-y[j]) == 1)) {
                    set.add_contact(i, j);
                }

            }
        }

        return set;

    } else if (dim == 3) {

        for (int i=0; i<n_res; i++) {
            for (int j=i+3; j<n_res; j++) {

                if (abs(j-i) % 2 == 0) continue;

                if ((abs(x[i]-x[j]) == 1 && abs(y[i]-y[j]) == 0 && abs(z[i]-z[j]) == 0) xor
                    (abs(x[i]-x[j]) == 0 && abs(y[i]-y[j]) == 1 && abs(z[i]-z[j]) == 0) xor
                    (abs(x[i]-x[j]) == 0 && abs(y[i]-y[j]) == 0 && abs(z[i]-z[j]) == 1)) {
                    set.add_contact(i, j);
                }

            }
        }

        return set;

    } else {

        std::cout << "Fatal error: dimension is not 2 or 3.\n";
        exit(1);

    }
}

