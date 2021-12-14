/*

    This file contains the basic
    objects employed by the rest 
    of the code. 

*/

#ifndef OBJECTS_H
#define OBJECTS_H

#include <list> 


class Sequence {

    private:

        int           len;
        int           dim;
        std::string   seq;

    public:

        Sequence(int len, int dim);

        void          set_seq(unsigned long long int k);
        std::string   get_seq();
        int           get_len();
        int           get_dim();
        void          reset();
        
};

class ContactSet {

    private:

        int               n_contacts;
        std::string       id;
        std::list <int>   first;
        std::list <int>   second;    

    public:

        ContactSet(std::string id);

        void              add_contact(int res1, int res2);

        int               get_n_contacts();
        std::string       get_id();
        std::list <int>   get_first_list();
        std::list <int>   get_second_list();

};

class Peptide {

    private:

        int               dim;
        int               n_res;
        int               res_placed;
        int*              x;
        int*              y;
        int*              z;

    public:

        Peptide(int n_res, int dim);
        ~Peptide();

        int               last_x();
        int               last_y();
        int               last_z();
        int               get_dim();
        int               get_seq_length();
        int               get_chain_length();

        bool              is_clashing(int _x, int _y, int _z=0);
        void              add_residue(int _x, int _y, int _z=0);
        ContactSet        get_contact_set(std::string seq);
        void              deallocate();
    
};

#endif
