#ifndef __BITSET_H__
#define __BITSET_H__

#include <cstdint>
#include <iostream>
#include <omp.h>

using namespace std;
class MyBitset{
    public:
    uint64_t* bits;
    int nbits;  // number of bits
    int nwords;
    int WORDSZ = 64;
    int num_1s = 0;

    MyBitset(int nbits){
        this->nbits = nbits;
        this->nwords = (nbits/WORDSZ) + 1; //idea: do ceil
        bits = new uint64_t[nwords];
        for (int i = 0; i < nwords; i++)
        {
            bits[i] = 0;
        }
    }

    //untested
    /*
    MyBitset(const MyBitset& other) {
        //     cout << "Copy constructor called" << endl;
        nbits = other.nbits;
        nwords = other.nwords;
        bits = new uint64_t[nwords];
        for (int i = 0; i < nwords; ++i) {
            bits[i] = other.bits[i];
        }
    }

    // Assignment operator
    MyBitset& operator=(const MyBitset& other) {
        if (this == &other) {
            return *this;
        }
        delete[] bits;

        nbits = other.nbits;
        nwords = other.nwords;

        bits = new uint64_t[nwords];
        for (int i = 0; i < nwords; ++i) {
            bits[i] = other.bits[i];
        }
        return *this;
    }
    // XOR operator
    MyBitset operator^(const MyBitset& other) const {
        if (nbits != other.nbits) {
            throw std::invalid_argument("Bitsets must be of the same size for XOR operation");
        }
        MyBitset result(nbits);
        for (int i = 0; i < nbits; ++i) {
            result.bits[i] = bits[i] ^ other.bits[i];
        }
        return result;
    }
    // MyBitset operator^(const MyBitset& b) { //wrong?
    //     MyBitset xor_bitset(this->nbits);

    //     //#pragma `omp parallel for
    //     for (int k = 0; k < this->nwords; k++) {
    //         xor_bitset.bits[k] = this->bits[k] ^ b.bits[k];
    //     }
    //     return xor_bitset;
    // }
    */

    int count_1s(){
        int sum = 0;
        //omp_set_num_threads(4);

        //#pragma `omp parallel for reduction(+:sum)
        for (int k = 0; k < nwords; k++) {
            sum += __builtin_popcountll ((uint64_t)bits[k]);
        }
        return sum;
    }

    vector<int> get_bits(){
        uint64_t bitset;
        vector<int> bitvec;
        for (int k = 0; k < nwords; k++) {
            bitset = bits[k];
            //std::cout<<"B"<<bitset<<std::endl;
            
            while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            int r = __builtin_ctzl(bitset);
            
            //std::cout<<(k * WORDSZ + r) << std::endl;
            
            
            bitvec.push_back(k * WORDSZ + r); //idea: reserve 1 counts
            //callback(k * WORDSZ + r);
            bitset ^= t;
            }
        }
        // for (int i = 0; i < nbits; i++)
        // {
        //     bitvec.push_back(get_bit(i));
        // }
        return bitvec;
    }

    void set_bit(int bit){
        bits[bit/WORDSZ] |= (uint64_t) 1 << (bit % WORDSZ);
    }

    void clear_bit(int bit){
        bits[bit/WORDSZ] &= ~((uint64_t) 1 << (bit % WORDSZ));
    }

    bool get_bit(int bit){
        return (bits[bit/WORDSZ] & ((uint64_t) 1 << (bit % WORDSZ))) != 0;
    }

    void print(){
        for (int i = 0; i < nbits; i++)
        {
            cout << get_bit(i);
        }
        cout << endl;
    }

    void print_pos(){
        for (int i = 0; i < nbits; i++)
        {
            if(get_bit(i)==1)
                cout << i << " ";
        }
        cout << endl;
    }

    ~MyBitset(){
        delete [] bits;
    }

};


#endif