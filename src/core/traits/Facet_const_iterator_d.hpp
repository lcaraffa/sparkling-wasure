/*
 * This file is part of the watertight surface reconstruction code https://github.com/lcaraffa/spark-ddt
 * Copyright (c) 2024 Caraffa Laurent, Mathieu Brédif.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef DDT_FACET_CONST_ITERATOR_D_HPP
#define DDT_FACET_CONST_ITERATOR_D_HPP

namespace ddt
{

template<typename TDS>
class Facet_const_iterator_d : public std::iterator<std::input_iterator_tag, typename TDS::Facet>
{
public:
    typedef typename TDS::Full_cell_const_handle     Cell_const_handle;
    typedef std::pair<Cell_const_handle, int> Facet;

    Facet_const_iterator_d()
        : tds_(nullptr), ft_(), cur_dim_(0)
    {
    }
    Facet_const_iterator_d(const TDS & tds)
        : tds_(&tds), ft_(tds.full_cells_begin(), 0), cur_dim_(tds.current_dimension())
    {
        assert( cur_dim_ > 0 );
        while( ! canonical() )
            raw_increment();
    }
    Facet_const_iterator_d(const TDS & tds, const Facet& ft)
        : tds_(&tds), ft_(ft), cur_dim_(tds.current_dimension())
    {
        assert( cur_dim_ > 0 );
        // do not enforce canonical here !
    }
    Facet_const_iterator_d(const Facet_const_iterator_d & fci)
        : tds_(fci.tds_), ft_(fci.ft_), cur_dim_(fci.cur_dim_)
    {
        assert( cur_dim_ > 0 );
        // do not enforce canonical here !
    }

    Facet_const_iterator_d & operator++()
    {
        increment();
        return (*this);
    }

    Facet_const_iterator_d operator++(int)
    {
        Facet_const_iterator_d tmp(*this);
        increment();
        return tmp;
    }

    bool operator==(const Facet_const_iterator_d & fi) const
    {
        if (tds_ == fi.tds_)
            return tds_ == nullptr ||
                   ((ft_.second == fi.ft_.second) && (ft_.first == fi.ft_.first));
        return (tds_ == nullptr && fi.tds_->full_cells_end() == fi.ft_.first ) ||
               (( fi.tds_ == nullptr && tds_->full_cells_end() == ft_.first ));
    }

    bool operator!=(const Facet_const_iterator_d & fi) const
    {
        return !(*this == fi);
    }

    const Facet& operator*() const
    {
        return ft_;
    }

    const Facet * operator->() const
    {
        return &ft_;
    }

    Facet_const_iterator_d& operator=(const Facet_const_iterator_d & fi)
    {
        tds_ = fi.tds_;
        ft_ = fi.ft_;
        cur_dim_ = fi.cur_dim_;
        return (*this);
    }

private:
    bool canonical()
    {
        if( tds_ == nullptr ) return true;
        if( tds_->full_cells_end() == ft_.first )
            return ( 0 == ft_.second );
        return ( ft_.first < ft_.first->neighbor(ft_.second) );
    }

    void raw_increment()
    {
        int i = ft_.second;
        if( i == cur_dim_ )
            ft_ = Facet(++ft_.first, 0);
        else
            ft_ = Facet(ft_.first, i + 1);
    }

    void increment()
    {
        do
        {
            raw_increment();
        }
        while( ! canonical() );
    }

    const TDS *tds_;
    Facet ft_;
    int cur_dim_;
};


}

#endif // DDT_FACET_CONST_ITERATOR_D_HPP
