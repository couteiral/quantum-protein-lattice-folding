/*

    Program to compute the contacts for
    all the possible conformations

*/

#include <cmath>
#include <iostream>
#include <fstream>

#include "objects.h"


Peptide check_seq(Sequence seq) {
    /*
        Reads a sequence and transforms it
        into a set of Cartesian coordinates.
    */

    Peptide pep(seq.get_len()+1, seq.get_dim());

    if (seq.get_dim() == 2) {

        int x, y;
        std::string sequence = seq.get_seq();

        for (int i=0; i<seq.get_len(); i++) {

            x = pep.last_x();
            y = pep.last_y();

            switch (sequence[i]) {

                case 'U': 
                    y++;
                    break;
                case 'D': 
                    y--;
                    break;
                case 'R': 
                    x++;
                    break;
                case 'L': 
                    x--;
                    break;

                default: 
                    std::cout << "Sequence contains wrong character " << sequence[i] << ".\n";
                    exit(1);

            }

            if (pep.is_clashing(x, y)) {
                return pep;
            } else {
                pep.add_residue(x, y);
            }

        }

        return pep;

    } else if (seq.get_dim() == 3) {

        int x, y, z;
        std::string sequence = seq.get_seq();

        for (int i=0; i<seq.get_len(); i++) {

            x = pep.last_x();
            y = pep.last_y();
            z = pep.last_z();

            switch (sequence[i]) {

                case 'U': 
                    y++;
                    break;
                case 'D': 
                    y--;
                    break;
                case 'R': 
                    x++;
                    break;
                case 'L': 
                    x--;
                    break;
                case 'I':
                    z++;
                    break;
                case 'O':
                    z--;
                    break;

                default: 
                    std::cout << "Sequence contains wrong character " << sequence[i] << ".\n";
                    exit(1);

            }

            if (pep.is_clashing(x, y, z)) {
                return pep;
            } else {
                pep.add_residue(x, y, z);
            }

        }

        return pep;

    } else {

        std::cout << "Fatal error: dimension is not 2 or 3.\n";
        exit(1);

    }

}

int comp_n_contacts(ContactSet cs1, ContactSet cs2) {

    if (cs1.get_n_contacts() > cs2.get_n_contacts()) {
        return true;
    } else {
        return false;
    }

}

int main(int argc, char *argv[]) {

    // Check the input arguments

    if (argc < 4) {
        std::cout << "Usage: \n";
        std::cout << "  cont_bf   seq_length  seq_dim  output_prefix \n";
        return 0;
    }

    int                    len = std::atoi(argv[1])-1;
    int                    dim = std::atoi(argv[2]);
    std::string            out(argv[3]);
    unsigned long long     max_combs;
    std::cout << max_combs << "\n";
    std::list <ContactSet> contacts;

    if (dim == 2) {
        max_combs = 1ULL << 2*len;
    } else if (dim == 3) {
        max_combs = (unsigned long long)(std::pow(6, len) + 0.5); // Only way, I know
    } else {
        std::cout << "ERROR: dimension can only be 2 or 3.\n";
    }

    std::cout << "Checking " << max_combs << " SAWs\n";

    // Iterate over all strings

    Sequence seq(len, dim);
    for (unsigned long long int k=0; k<max_combs; k++) {

        seq.set_seq(k);
        //std::cout << "   checking " << seq.get_seq() << "...\n";

        Peptide pep = check_seq(seq);

        if (pep.get_chain_length() == pep.get_seq_length()) {
            ContactSet set = pep.get_contact_set(seq.get_seq());
            if (set.get_n_contacts() != 0) contacts.push_back(set);
        }

        seq.reset();

    }

    // Print everything

    std::ofstream seq_file;
    std::ofstream contact_file;
    std::ofstream meta_file;
    
    seq_file.open(out + ".seq");
    contact_file.open(out + ".contacts");

    contacts.sort(comp_n_contacts); // Sort in descending order of contacts

    meta_file.open(out + ".meta");
    meta_file << contacts.size();
    meta_file.close();

    for (auto it : contacts) {

        seq_file << it.get_id() << "\n";

        for (int el1 : it.get_first_list()) {
            contact_file << el1 << " ";
        }
        for (int el2 : it.get_second_list()) {
            contact_file << el2 << " ";
        }
        contact_file << "\n";

    }

    seq_file.close();
    contact_file.close();

    return 0;
        
}
