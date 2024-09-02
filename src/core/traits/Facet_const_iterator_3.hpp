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
#ifndef DDT_FACET_CONST_ITERATOR_3_HPP
#define DDT_FACET_CONST_ITERATOR_3_HPP

namespace ddt
{

template <typename TDS> using Cell_const_iterator_3 = typename TDS::Cell_iterator;
template <typename TDS> using Facet_3 = std::pair<Cell_const_iterator_3<TDS>, int>;

template<typename TDS>
class Facet_const_iterator_3 : public std::iterator<std::input_iterator_tag, Facet_3<TDS>>
{
public:
    typedef Cell_const_iterator_3<TDS> Cell_const_iterator;
    typedef Facet_3<TDS>               Facet;

    Facet_const_iterator_3()
        : tds_(nullptr), ft_()
    {
    }
    Facet_const_iterator_3(const TDS & tds)
        : tds_(&tds), ft_(tds.cells_begin(), 0)
    {
        while( ! canonical() ) raw_increment();
    }
    Facet_const_iterator_3(const TDS & tds, const Facet& ft)
        : tds_(&tds), ft_(ft)
    {
        // do not enforce canonical here !
    }
    Facet_const_iterator_3(const Facet_const_iterator_3 & fci)
        : tds_(fci.tds_), ft_(fci.ft_)
    {
        // do not enforce canonical here !
    }

    Facet_const_iterator_3 & operator++()
    {
        increment();
        return (*this);
    }

    Facet_const_iterator_3 operator++(int)
    {
        Facet_const_iterator_3 tmp(*this);
        increment();
        return tmp;
    }

    bool operator==(const Facet_const_iterator_3 & fi) const
    {
        if (tds_ == fi.tds_)
            return tds_ == nullptr ||
                   ((ft_.second == fi.ft_.second) && (ft_.first == fi.ft_.first));
        return (tds_ == nullptr && fi.tds_->cells_end() == fi.ft_.first ) ||
               (( fi.tds_ == nullptr && tds_->cells_end() == ft_.first ));
    }

    bool operator!=(const Facet_const_iterator_3 & fi) const
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

    Facet_const_iterator_3& operator=(const Facet_const_iterator_3 & fi)
    {
        tds_ = fi.tds_;
        ft_ = fi.ft_;
        return (*this);
    }

private:
    bool canonical()
    {
        if( tds_ == nullptr ) return true;
        if( tds_->cells_end() == ft_.first )
            return ( 0 == ft_.second );
        return ( ft_.first < ft_.first->neighbor(ft_.second) );
    }

    void raw_increment()
    {
        int i = ft_.second;
        if( i == 3 )
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
};


}

#endif // DDT_FACET_CONST_ITERATOR_3_HPP
