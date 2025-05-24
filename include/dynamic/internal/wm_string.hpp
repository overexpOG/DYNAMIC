// Copyright (c) 2019, Erik Garrison.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * wm_string.hpp
 *
 *  Created on: May 30, 2020
 *      Author: MitI_7
 *      Author: Erik Garrison
 *
 *  Dynamic string supporting rank, select, access, insert.
 *  Backed by a wavelet matrix based on the templated dynamic bitvector.
 *
 */

#ifndef INCLUDE_INTERNAL_WM_STRING_HPP_
#define INCLUDE_INTERNAL_WM_STRING_HPP_

#include "dynamic/internal/includes.hpp"

namespace dyn
{

    template <class dynamic_bitvector_t>

    class wm_string
    {

    public:
        std::vector<dynamic_bitvector_t> bit_arrays; // Array of bitvectors
        std::vector<ulint> begin_one;                // Initial position of 1 for each bit
                                                     // Amount of zeroes in each level

        ulint n = 0;         // size of the matrix
        ulint sigma = 0;     // max char
        ulint bit_width = 0; // number of bits necessary to represent a character

    public:
        wm_string(void) {} // default constructor, for loading from file

        // max_element: the biggest number that enters
        wm_string(ulint sigma) : n(0), sigma(sigma + 1)
        {
            this->bit_width = this->get_num_of_bit(sigma);
            if (bit_width == 0)
            {
                bit_width = 1;
            }
            this->begin_one.resize(bit_width);

            for (ulint i = 0; i < bit_width; ++i)
            {
                dynamic_bitvector_t sv;
                bit_arrays.push_back(sv);
            }
        }

        wm_string(ulint num_of_alphabet, const std::vector<ulint> &array) : n(0), sigma(num_of_alphabet + 1)
        {
            this->bit_width = this->get_num_of_bit(num_of_alphabet);
            if (bit_width == 0)
            {
                bit_width = 1;
            }
            this->begin_one.resize(bit_width);

            if (array.empty())
            {
                for (ulint i = 0; i < bit_width; ++i)
                {
                    dynamic_bitvector_t sv;
                    bit_arrays.push_back(sv);
                }
                return;
            }

            n = array.size();

            std::vector<ulint> v(array), b(array.size(), 0);

            for (ulint i = 0; i < bit_width; ++i)
            {

                std::vector<ulint> temp;
                // put 0 in temp
                for (ulint j = 0; j < v.size(); ++j)
                {
                    ulint c = v.at(j);
                    ulint bit = (c >> (bit_width - i - 1)) & 1; // i-th bit from the top
                    if (bit == 0)
                    {
                        temp.push_back(c);
                        b[j] = 0;
                    }
                }

                this->begin_one.at(i) = temp.size();

                // put 1 in temp
                for (ulint j = 0; j < v.size(); ++j)
                {
                    ulint c = v.at(j);
                    ulint bit = (c >> (bit_width - i - 1)) & 1; // 　i-th bit from the top
                    if (bit == 1)
                    {
                        temp.push_back(c);
                        b[j] = 1;
                    }
                }

                dynamic_bitvector_t dbv;
                for (auto &i : b)
                {
                    dbv.push_back(i);
                }
                bit_arrays.emplace_back(dbv);
                v = temp;
            }
        }

        ulint serialize(ostream &out)
        {
            ulint w_bytes = 0;
            out.write((char *)&n, sizeof(n));
            w_bytes += sizeof(n);
            out.write((char *)&sigma, sizeof(sigma));
            w_bytes += sizeof(sigma);
            out.write((char *)&bit_width, sizeof(bit_width));
            w_bytes += sizeof(bit_width);
            out.write((char *)begin_one.data(), sizeof(ulint) * bit_width);
            w_bytes += sizeof(ulint) * bit_width;
            for (auto &bv : bit_arrays)
            {
                w_bytes += bv.serialize(out);
            }
            return w_bytes;
        }

        void load(istream &in)
        {
            in.read((char *)&n, sizeof(n));
            in.read((char *)&sigma, sizeof(sigma));
            in.read((char *)&bit_width, sizeof(bit_width));
            begin_one.resize(bit_width);
            in.read((char *)begin_one.data(), sizeof(ulint) * bit_width);
            bit_arrays.resize(bit_width);
            for (ulint i = 0; i < bit_width; ++i)
            {
                bit_arrays[i].load(in);
            }
        }

        // v[pos]
        ulint at(ulint pos)
        {
            assert(pos < this->n);

            ulint c = 0;
            for (ulint i = 0; i < bit_arrays.size(); ++i)
            {
                dynamic_bitvector_t &level = bit_arrays.at(i);
                ulint bit = level.at(pos); // i-th bit of the original value
                c = (c <<= 1) | bit;
                pos = level.rank(pos, bit);
                if (bit)
                {
                    pos += this->begin_one.at(i);
                }
            }
            return c;
        }

        // high-level access
        ulint operator[](ulint i) { return this->at(i); }

        void increment_alphabet()
        {
            bit_width = this->get_num_of_bit(++this->sigma - 1);
            if (bit_arrays.size() < bit_width)
            {
                dynamic_bitvector_t dbv;
                for (uint64_t i = 0; i < this->n; i++)
                {
                    dbv.push_back(0);
                }
                bit_arrays.emplace(bit_arrays.begin(), dbv);
                this->begin_one.emplace(this->begin_one.begin(), this->n);
            }
        }

        // number of c's in v[0..pos)
        ulint rank(ulint pos, ulint c)
        {
            assert(pos <= n);
            if (c >= sigma)
            {
                return 0;
            }

            ulint left = 0, right = pos;
            for (ulint i = 0; i < bit_width; ++i)
            {
                const ulint bit = (c >> (bit_width - i - 1)) & 1; // i-th bit from the top
                dynamic_bitvector_t &level = bit_arrays.at(i);
                left = level.rank(left, bit);   // number of numbers that are equal to the i-th bit of c
                right = level.rank(right, bit); // number of numbers that are equal to the i-th bit of c
                if (bit)
                {
                    left += this->begin_one.at(i);
                    right += this->begin_one.at(i);
                }
            }

            return right - left;
        }

        // Returns the i-th position of c+1
        // The rank is 1-origin
        ulint select(ulint rank, ulint c)
        {
            assert(rank > 0);
            assert(c < sigma);
            --rank;

            ulint left = 0;
            for (ulint i = 0; i < bit_width; ++i)
            {
                const ulint bit = (c >> (bit_width - i - 1)) & 1; // i-th bit from the top
                left = bit_arrays.at(i).rank(left, bit);          // number of numbers that are equal to the i-th bit of c
                if (bit)
                {
                    left += this->begin_one.at(i);
                }
            }

            ulint index = left + rank;
            for (ulint i = 0; i < bit_arrays.size(); ++i)
            {
                ulint bit = ((c >> i) & 1); // i-th bit from the bottom
                if (bit == 1)
                {
                    index -= this->begin_one.at(bit_width - i - 1);
                }
                index = this->bit_arrays.at(bit_width - i - 1).select(index, bit);
            }
            return index + 1;
        }

        // Insert c in pos
        void insert(ulint pos, ulint c)
        {
            assert(pos <= this->n);

            for (ulint i = 0; i < bit_arrays.size(); ++i)
            {
                const ulint bit = (c >> (bit_width - i - 1)) & 1; // i-th bit from the top
                bit_arrays.at(i).insert(pos, bit);
                pos = bit_arrays.at(i).rank(pos, bit);
                if (bit)
                {
                    pos += this->begin_one.at(i);
                }
                else
                {
                    this->begin_one.at(i)++;
                }
            }

            this->n++;
        }

        void push_front(ulint c)
        {
            this->insert(0, c);
        }

        // add c ath the end
        void push_back(ulint c)
        {
            this->insert(this->n, c);
        }

        // remove pos
        void remove(ulint pos)
        {
            assert(pos < this->n);
            if (pos >= this->n)
            {
                throw "Segmentation fault";
            }

            for (ulint i = 0; i < bit_arrays.size(); ++i)
            {
                ulint bit = bit_arrays.at(i).at(pos); // the i-th bit of the original number

                auto next_pos = bit_arrays.at(i).rank(pos, bit);
                bit_arrays.at(i).remove(pos);

                if (bit)
                {
                    next_pos += this->begin_one.at(i);
                }
                else
                {
                    this->begin_one.at(i)--;
                }
                pos = next_pos;
            }
            this->n--;
        }

        uint64_t remove_and_return(uint64_t pos)
        {
            assert(pos < this->n);
            if (pos >= this->n)
            {
                throw "Segmentation fault";
            }

            uint64_t res = 0;

            for (ulint i = 0; i < bit_arrays.size(); ++i)
            {
                ulint bit = bit_arrays.at(i).at(pos); // the i-th bit of the original number
                res <<= 1;
                res |= bit;

                auto next_pos = bit_arrays.at(i).rank(pos, bit);
                bit_arrays.at(i).remove(pos);

                if (bit)
                {
                    next_pos += this->begin_one.at(i);
                }
                else
                {
                    this->begin_one.at(i)--;
                }
                pos = next_pos;
            }
            this->n--;
            return res;
        }

        // update c in pos
        void update(ulint pos, ulint c)
        {
            assert(pos < this->n);
            this->remove(pos);
            this->insert(pos, c);
        }

        ulint size(void) const
        {
            return this->n;
        }

        ulint bit_size(void) const
        {
            ulint n_bits = 0;
            for (auto &ba : bit_arrays)
            {
                n_bits += ba.bit_size();
            }
            n_bits += sizeof(ulint) * begin_one.size();
            return n_bits;
        }

        // Other operations are the same as the normal wavelet matrix

        /*** Select next
         * Implemented by Yuval Linker
         * Find the position of the next c outside of the range [0,pos] 
         * 
         * @arg i
         * @arg c
         * @arg n_elems
         * @returns pair<>
         * */
        std::pair<ulint, ulint> select_next(ulint pos, ulint c, ulint n_elems)
        {
            assert(pos <= size());
            uint64_t mask = 1ULL << (bit_arrays.size() - 1);
            // Interval [b...b+r]
            ulint b = 0;
            ulint r = pos; // size of interval
            ulint rank_b = 0;

            // Vectors
            std::vector<ulint> path_beginnings(bit_arrays.size() + 1);
            std::vector<ulint> path_ranks(bit_arrays.size() + 1);

            for (ulint k = 0; k < bit_arrays.size() and pos; k++)
            {
                dynamic_bitvector_t &level = bit_arrays.at(k);
                rank_b = level.rank(b);                  // ones in [0...b)
                ulint ones = level.rank(b + r) - rank_b; // ones in [b...r)
                if (c & mask)
                {
                    r = ones;
                    b = begin_one.at(k) + rank_b;
                }
                else
                {
                    r = r - ones;
                    b -= rank_b;
                }
                mask >>= 1;
                path_beginnings[k + 1] = b;
                path_ranks[k] = rank_b;
            }

            mask = 1ULL;
            pos = r + 1;
            if (pos > n_elems)
                return std::pair<ulint, ulint>(0, 0);
            uint64_t levels_amount = bit_arrays.size();

            for (ulint k = 1; k <= levels_amount; k++)
            {
                b = path_beginnings[levels_amount - k];
                rank_b = path_ranks[levels_amount - k];
                dynamic_bitvector_t &level = bit_arrays.at(levels_amount - k);
                if (c & mask)
                { // right child => search i'th one
                    pos = level.select1(rank_b + pos - 1) - b + 1;
                }
                else
                { // left child => search i'th zero
                    pos = level.select0(b - rank_b + pos - 1) - b + 1;
                }
                mask <<= 1;
            }
            return std::pair<ulint, ulint>(pos - 1, r);
        }

        ulint range_minimum_query(ulint i, ulint j)
        {
            return _range_minimum_query(i, j, 0, 0, 0);
        }

        ulint _range_minimum_query(ulint i, ulint j, ulint depth, ulint b, ulint res)
        {
            if (depth == bit_arrays.size())
            {
                return res;
            }
            else
            {
                dynamic_bitvector_t &level = bit_arrays.at(depth);
                ulint rank_0_b = level.rank(b);
                ulint rank_b_i = level.rank(b + i) - rank_0_b;
                ulint rank_b_j = level.rank(b + j + 1) - rank_0_b;

                ulint i_l = i - rank_b_i;
                ulint j_l = j - rank_b_j;
                ulint i_r = i - i_l;
                ulint j_r = j - 1 - j_l;
                ulint n_l = j_l - i_l + 1;

                res <<= 1;
                if (n_l == 0)
                {
                    b = begin_one.at(depth) + rank_0_b;
                    res |= 1;
                    return _range_minimum_query(i_r, j_r, depth + 1, b, res);
                }
                else
                {
                    b -= rank_0_b;
                    return _range_minimum_query(i_l, j_l, depth + 1, b, res);
                }
            }
        }

        // TODO: Document
        /*** Recursively finds the value in range that is >= x
         * \param
         *
         * */
        ulint range_next_value(ulint x, ulint i, ulint j) 
        {
            // c is greater than any symbol in the wm
            if (((1ULL) << bit_arrays.size()) <= x)
                return 0;

            return _range_next_value(x, i, j, 0, 0, 0);
        }

        ulint _range_next_value(ulint x, ulint i, ulint j, ulint depth, ulint b, ulint res) 
        {
            if (b + i > b + j)
            {
                return 0;
            }
            else
            {
                if (depth == bit_arrays.size())
                    return res;
                else
                {
                    dynamic_bitvector_t &level = bit_arrays.at(depth);
                    ulint rank_0_b = level.rank(b);
                    ulint rank_b_i = level.rank(b + i) - rank_0_b;
                    ulint rank_b_j = level.rank(b + j + 1) - rank_0_b;

                    ulint i_l = i - rank_b_i; // i in the left child
                    ulint j_l = j - rank_b_j; // j in the left child
                    ulint i_r = i - i_l;      // i in the right child
                    ulint j_r = j - 1 - j_l;  // j in the right child

                    ulint mask = (1ULL) << (bit_arrays.size() - 1 - depth);
                    res <<= 1;
                    if (x & mask)
                    { // recurse on the right child
                        b = begin_one.at(depth) + rank_0_b;
                        res |= 1;
                        return _range_next_value(x, i_r, j_r, depth + 1, b, res);
                    }
                    else
                    { // recurse on the left child
                        b -= rank_0_b;
                        ulint y = _range_next_value(x, i_l, j_l, depth + 1, b, res);

                        if (y != 0)
                            return y;
                        else
                        {
                            // It didnt find the next value in the left child so find the min in the right child
                            b = begin_one.at(depth) + rank_0_b;
                            res |= 1;
                            return _range_next_value_min(i_r, j_r, depth + 1, b, res);
                        }
                    }
                }
            }
        }

        ulint _range_next_value_min(ulint i, ulint j, ulint depth, ulint b, ulint res)
        {
            if (b + i > b + j)
                return 0;

            if (depth == bit_arrays.size())
                return res;
            else
            {
                dynamic_bitvector_t &level = bit_arrays.at(depth);
                ulint rank_0_b = level.rank(b);
                ulint rank_b_i = level.rank(b + i) - rank_0_b;
                ulint rank_b_j = level.rank(b + j + 1) - rank_0_b;

                ulint i_l = i - rank_b_i;
                ulint j_l = j - rank_b_j;
                ulint i_r = i - i_l;
                ulint j_r = j - 1 - j_l;
                ulint n_l = j_l - i_l + 1; // number of elements in the interval on the left child

                res <<= 1;
                if (n_l == 0)
                { // recurse on the right child
                    b = begin_one.at(depth) + rank_0_b;
                    res |= 1;
                    return _range_next_value_min(i_r, j_r, depth + 1, b, res);
                }
                else
                { // recurse on the left child
                    b -= rank_0_b;
                    return _range_next_value_min(i_l, j_l, depth + 1, b, res);
                }
            }
        }

        std::pair<range_t, range_t> child_ranges(ulint depth, range_t r)
        {
            dynamic_bitvector_t &level = bit_arrays.at(depth);
            ulint left_of_range_ones = level.rank(r.first);
            ulint in_range_ones = level.rank(r.second + 1);
            ulint right_size = in_range_ones - left_of_range_ones;
            ulint left_size = (r.second - r.first + 1) - right_size;

            ulint left_offset = r.first - left_of_range_ones;
            ulint right_offset = left_of_range_ones + begin_one.at(depth);

            return {
                {left_offset, left_offset + left_size - 1},
                {right_offset, right_offset + right_size - 1}};
        }

        std::vector<ulint>
        all_values_in_range(ulint lb, ulint rb)
        {
            std::vector<ulint> res_vec;
            if (lb <= rb)
            {
                _all_values_in_range(0, 0, {lb, rb}, res_vec);
            }
            return res_vec;
        }

        void _all_values_in_range(ulint acc_sym, ulint depth, range_t r, std::vector<ulint> &res_vec)
        {
            using std::get;
            if (get<0>(r) > get<1>(r))
                return;

            ulint level_sym = acc_sym << 1;

            if (depth == bit_arrays.size() - 1)
            {
                ulint l_bound = get<0>(r);
                ulint r_bound = get<1>(r);
                dynamic_bitvector_t &level = bit_arrays.at(depth);
                ulint rank_outside_range = level.rank(l_bound);
                ulint rank_inside_range = level.rank(r_bound + 1);
                if ((r_bound - l_bound + 1 - (rank_inside_range - rank_outside_range)) > 0)
                {
                    res_vec.emplace_back(level_sym);
                }
                if (rank_inside_range - rank_outside_range > 0)
                {
                    level_sym |= 1;
                    res_vec.emplace_back(level_sym);
                }
                return;
            }

            // Transform range to left and right child

            std::pair<range_t, range_t> next_ranges = child_ranges(depth, r);

            if (!dyn::empty(get<0>(next_ranges)))
            { // Left child
                _all_values_in_range(level_sym, depth + 1, get<0>(next_ranges), res_vec);
            }
            if (!dyn::empty(get<1>(next_ranges)))
            { // right child
                _all_values_in_range(level_sym | 1, depth + 1, get<1>(next_ranges), res_vec);
            }
        }

        std::pair<ulint, ulint> inverse_select(ulint i)
        {
            assert(i < size());
            ulint c = 0;
            ulint b = 0; // start position of the interval
            ulint mask = (1ULL) << (bit_arrays.size() - 1);
            for (ulint k = 0; k < bit_arrays.size(); ++k)
            {
                dynamic_bitvector_t &level = bit_arrays.at(k);
                ulint ones_left_of_interval = level.rank(b);
                ulint ones_in_interval = level.rank(b + i) - ones_left_of_interval;
                c <<= 1;
                if (level[b + i])
                {
                    i = ones_in_interval;
                    b = begin_one.at(k) + ones_left_of_interval;
                    c |= 1;
                }
                else
                {
                    i -= ones_in_interval;
                    b -= ones_left_of_interval;
                }
                mask >>= 1;
            }
            return {i, c};
        }

    private:
        ulint get_num_of_bit(ulint x)
        {
            if (x == 0)
                return 0;

            return 64 - __builtin_clzll(x);
        }
    };

}

#endif /* INCLUDE_INTERNAL_WM_STRING_HPP_ */
