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
#ifndef DDT_CELL_CONST_ITERATOR_HPP
#define DDT_CELL_CONST_ITERATOR_HPP

#include "ddt_exeptions.hpp"

namespace ddt
{

template<typename DDT>
class Cell_const_iterator : public std::iterator<std::input_iterator_tag, Cell_const_iterator<DDT> >
{
public:
    typedef typename DDT::Traits                    Traits;
    typedef typename Traits::Flag_C                    Flag_C;
    typedef typename Traits::Data_C                    Data_C;
    typedef typename DDT::Id                        Id;
    typedef typename DDT::Tile_cell_const_handle    Tile_cell_const_handle;
    typedef typename DDT::Tile_cell_const_iterator  Tile_cell_const_iterator;
    typedef typename DDT::Tile_vertex_const_handle  Tile_vertex_const_handle;
    typedef typename DDT::Tile_const_iterator       Tile_const_iterator;
    typedef typename DDT::Vertex_const_iterator     Vertex_const_iterator;
    typedef typename DDT::Facet_const_iterator      Facet_const_iterator;

private:
    Tile_const_iterator begin_;
    Tile_const_iterator end_;
    Tile_const_iterator tile_;
    Tile_cell_const_iterator cell_;

public:
    Cell_const_iterator(Tile_const_iterator begin, Tile_const_iterator end)
        : begin_(begin), end_(end), tile_(begin), cell_()
    {
        if(tile_ != end_)
        {
            cell_ = tile_->cells_begin();
            advance_to_main();
        }
        assert(is_valid());
    }

    Cell_const_iterator(Tile_const_iterator begin, Tile_const_iterator end, Tile_const_iterator tile)
        : begin_(begin), end_(end), tile_(tile), cell_()
    {
        if(tile_ != end_)
        {
            cell_ = tile_->cells_begin();
            advance_to_main();
        }
        assert(is_valid());
    }

    Cell_const_iterator(Tile_const_iterator begin, Tile_const_iterator end, Tile_const_iterator tile, Tile_cell_const_iterator cell)
        : begin_(begin), end_(end), tile_(tile), cell_(cell)
    {
        // do not enforce main here !
        assert(is_valid());
    }

    Cell_const_iterator(const Cell_const_iterator& c)
        : begin_(c.begin_), end_(c.end_), tile_(c.tile_), cell_(c.cell_)
    {
        // do not enforce main here !
        assert(is_valid());
    }

    Cell_const_iterator& advance_to_main()
    {
        while(tile_ != end_)
        {
            if(cell_ == tile_->cells_end())
            {
                if (++tile_ != end_) cell_ = tile_->cells_begin();
            }
            else if(tile_->cell_is_main(cell_))
            {
                break;
            }
            else
            {
                ++cell_;
            }
        }
        return *this;
    }

    bool operator<(const Cell_const_iterator& c) const
    {
        if (c.tile_ == c.end_) return tile_ != end_;
        if (tile_ == end_) return false;
        return  tile_->id() < c.tile_->id() || (tile_->id() == c.tile_->id() && cell_ < c.cell_);
    }

    Cell_const_iterator& operator++()
    {
        assert(tile_ != end_);
        ++cell_;
        return advance_to_main();
    }

    Cell_const_iterator operator++(int)
    {
        Cell_const_iterator tmp(*this);
        operator++();
        return tmp;
    }

    Cell_const_iterator& operator+=(int n)
    {
        assert(tile_ != end_);
        for(auto cit = tile_->cells_begin(); cit != cell_; ++cit)
            if(tile_->cell_is_main(cit))
                ++n;
        int num_main_cells = tile_->number_of_main_cells();
        while(n >= num_main_cells)
        {
            n -= num_main_cells;
            ++tile_;
            num_main_cells = tile_->number_of_main_cells();
        }
        cell_ = tile_->cells_begin();
        advance_to_main();
        for(; n>0 ; --n)
            ++(*this);
        assert(is_valid());
        return *this;
    }

    bool operator==(const Cell_const_iterator& rhs) const
    {
        return begin_ == rhs.begin_
               && tile_ == rhs.tile_
               && end_ == rhs.end_
               && (tile_ == end_ || cell_==rhs.cell_);
    }

    bool operator!=(const Cell_const_iterator& rhs) const { return !(*this == rhs); }
    Cell_const_iterator& operator*() { return *this; }
    Cell_const_iterator* operator->() { return this; }

    Facet_const_iterator facet(int i) const
    {
        assert(tile_ != end_);
        return Facet_const_iterator(begin_, end_, tile_, tile_->facet(cell_, i));
    }

    Cell_const_iterator neighbor(int i) const
    {
        assert(tile_ != end_);
        Tile_cell_const_iterator c = cell_->neighbor(i);
        if(!tile_->cell_is_foreign(c))
            return Cell_const_iterator(begin_, end_, tile_, c);
        // there is no representative of the neighbor in tile_
        return facet(i)->main()->neighbor()->full_cell();
    }




    int mirror_index(int i) const
    {
        assert(tile_ != end_);
        Tile_cell_const_iterator c = cell_->neighbor(i);
        if(!tile_->cell_is_foreign(c))
            return tile_->mirror_index(cell_, i);
        return facet(i)->main()->mirror_index();
    }

    Id main_id() const
    {
        assert(tile_ != end_);
        return tile_->cell_main_id(cell_);
    }


    Id tile_id() const
    {
        assert(tile_ != end_);
        return tile_->id();
    }


    Id gid() const
    {
        assert(tile_ != end_);
        return tile_->gid(cell_);
    }


    Id lid() const
    {
        assert(tile_ != end_);
        return tile_->lid(cell_);
    }



    Vertex_const_iterator vertex (const int i) const
    {
        assert(tile_ != end_);
        return Vertex_const_iterator(begin_, end_, tile_, tile_->vertex(cell_, i));
    }


    const  Data_C & cell_data() const
    {
        return tile_->datac(cell_);
    }



    Flag_C & flag()
    {
        return tile_->flagc(cell_);
    }



    Cell_const_iterator main() const
    {
        assert(tile_ != end_);
        Id id = main_id();
        if (id == tile_->id()) return *this; // <=> is_main
        for (Tile_const_iterator tile = begin_; tile != end_; ++tile)
            if (id == tile->id())
            {
                Tile_cell_const_iterator c = tile->locate_cell(*tile_, cell_);
                if (c==tile->cells_end()) return Cell_const_iterator(begin_, end_, end_);
                return Cell_const_iterator(begin_, end_, tile, c);
            }
        throw DDT_exeption("ex_main_cell_not_found");
        assert(false);
        // no tile has the main id in tiles
        return Cell_const_iterator(begin_, end_, end_);
    }



    Tile_const_iterator    tile()      const { return tile_; }
    Tile_cell_const_handle full_cell() const { return cell_; }

    bool is_local()    const { return tile_->cell_is_local(cell_); }
    int local_score()    const { return tile_->cell_local_score(cell_); }
    bool is_mixed()    const { return tile_->cell_is_mixed(cell_); }
    bool is_foreign()  const { return tile_->cell_is_foreign(cell_); }
    bool is_main()     const { return tile_->cell_is_main(cell_); }
    bool has_id(Id ii)     const { return tile_->has_id(cell_,ii); }
    bool is_infinite() const { return tile_->cell_is_infinite(cell_); }
    bool is_inside() const { return tile_->cell_is_inside(cell_); }

    bool circumcircle() const { return tile_->cell_is_inside(cell_); }
    std::vector<double> barycenter() const { return tile_->get_cell_barycenter(cell_); }

    void get_list_vertices(std::list<Tile_vertex_const_handle> & ll) const { return tile_->get_list_vertices(cell_,ll);}

    bool is_valid()    const
    {
        return tile_ == end_ || cell_ != tile_->cells_end();
    }
};

}

#endif // DDT_CELL_CONST_ITERATOR_HPP







