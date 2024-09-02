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
#ifndef DDT_FACET_CONST_ITERATOR_HPP
#define DDT_FACET_CONST_ITERATOR_HPP

#include "ddt_exeptions.hpp"

namespace ddt
{

template<typename DDT>
class Facet_const_iterator : public std::iterator<std::input_iterator_tag, Facet_const_iterator<DDT> >
{
public:
    typedef typename DDT::Traits                    Traits;
    typedef typename DDT::Id                        Id;
    typedef typename DDT::Tile_facet_const_iterator Tile_facet_const_iterator;
    typedef typename DDT::Tile_cell_const_iterator  Tile_cell_const_iterator;
    typedef typename DDT::Tile_const_iterator       Tile_const_iterator;
    typedef typename DDT::Cell_const_iterator       Cell_const_iterator;

private:
    Tile_const_iterator begin_;
    Tile_const_iterator end_;
    Tile_const_iterator tile_;
    Tile_facet_const_iterator facet_;

public:
    Facet_const_iterator(Tile_const_iterator begin, Tile_const_iterator end)
        : begin_(begin), end_(end), tile_(begin), facet_()
    {
        if(tile_ != end_)
        {
            facet_ = tile_->facets_begin();
            advance_to_main();
        }
        assert(is_valid());
    }

    Facet_const_iterator(Tile_const_iterator begin, Tile_const_iterator end, Tile_const_iterator tile)
        : begin_(begin), end_(end), tile_(tile), facet_()
    {
        if(tile_ != end_)
        {
            facet_ = tile_->facets_begin();
            advance_to_main();
        }
        assert(is_valid());
    }

    Facet_const_iterator(Tile_const_iterator begin, Tile_const_iterator end, Tile_const_iterator tile, Tile_facet_const_iterator facet)
        : begin_(begin), end_(end), tile_(tile), facet_(facet)
    {
        // do not enforce main here !
        assert(is_valid());
    }

    Facet_const_iterator(const Facet_const_iterator& c)
        : begin_(c.begin_), end_(c.end_), tile_(c.tile_), facet_(c.facet_)
    {
        // do not enforce main here !
        assert(is_valid());
    }
    Facet_const_iterator& advance_to_main()
    {
        while(tile_ != end_)
        {
            if(facet_ == tile_->facets_end())
            {
                if (++tile_ != end_) facet_ = tile_->facets_begin();
            }
            else if(tile_->facet_is_main(facet_))
            {
                break;
            }
            else
            {
                ++facet_;
            }
        }
        return *this;
    }

    Facet_const_iterator& operator++()
    {
        assert(tile_ != end_);
        ++facet_;
        return advance_to_main();
    }

    Facet_const_iterator operator++(int)
    {
        Facet_const_iterator tmp(*this);
        operator++();
        return tmp;
    }

    bool operator==(const Facet_const_iterator& rhs) const
    {
        return begin_ == rhs.begin_
               && tile_ == rhs.tile_
               && end_ == rhs.end_
               && (tile_ == end_ || facet_ == rhs.facet_);
    }

    bool operator!=(const Facet_const_iterator& rhs) const { return !(*this == rhs); }
    Facet_const_iterator& operator*() { return *this; }
    Facet_const_iterator* operator->() { return this; }

    Id main_id() const
    {
        assert(tile_ != end_);
        Id id = -1;
        int cid = index_of_covertex();
        auto c = tile_->full_cell(facet_);
        int D = tile_->current_dimension();
        for(int i=0; i<=D; ++i)
        {
            if (i == cid) continue;
            auto v = tile_->vertex(c, i);
            if (tile_->vertex_is_infinite(v)) continue;
            Id vid = tile_->id(v);
            if (id == -1 || vid < id) id = vid;
        }
        return id;
    }



    Facet_const_iterator main() const
    {
        assert(tile_ != end_);
        Id id = main_id();
        if (id == tile_->id()) // <=> is_main
            return *this;
        for (Tile_const_iterator t = begin_; t != end_; ++t)
            if (id == t->id())
                return Facet_const_iterator(begin_, end_, t, t->locate_facet(*tile_, facet_));
        std::cerr << "catch neighbor exeption" << std::endl;
        throw DDT_exeption("ex_main_facet_not_found");
        assert(false);
        return Facet_const_iterator(begin_, end_, end_);
    }

    Facet_const_iterator neighbor() const
    {
        assert(tile_ != end_);
        int i = tile_->index_of_covertex(facet_);
        Tile_cell_const_iterator c = tile_->full_cell(facet_);
        Tile_cell_const_iterator n = tile_->neighbor(c, i);
        if (!tile_->cell_is_foreign(n))
            return Facet_const_iterator(begin_, end_, tile_, tile_->facet(n, tile_->mirror_index(c,i)));
        return main()->neighbor();
    }

    inline int mirror_index() const
    {
        return main()->neighbor()->index_of_covertex();
    }

    inline Cell_const_iterator full_cell() const
    {
        assert(tile_ != end_);
        return Cell_const_iterator(begin_, end_, tile_, tile_->full_cell(facet_));
    }

    inline int index_of_covertex() const
    {
        assert(tile_ != end_);
        return tile_->index_of_covertex(facet_);
    }

    Tile_const_iterator       tile()  const { return tile_;  }
    Tile_facet_const_iterator facet() const { return facet_; }

    bool is_local()    const { return tile_->facet_is_local(facet_); }
    int local_score() const { return tile_->facet_local_score(facet_); }
    bool is_mixed()    const { return tile_->facet_is_mixed(facet_); }
    bool facet_has_id(Id ii)    const { return tile_->facet_has_id(facet_,ii); }
    bool is_foreign()  const { return tile_->facet_is_foreign(facet_); }
    bool is_main()     const { return tile_->facet_is_main(facet_); }
    bool is_infinite() const { return tile_->facet_is_infinite(facet_); }

    bool is_valid()    const
    {
        return tile_ == end_ || facet_ != tile_->facets_end();
    }
};

}

#endif // DDT_FACET_CONST_ITERATOR_HPP
