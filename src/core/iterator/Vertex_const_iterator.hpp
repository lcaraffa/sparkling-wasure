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
#ifndef DDT_VERTEX_CONST_ITERATOR_HPP
#define DDT_VERTEX_CONST_ITERATOR_HPP

namespace ddt
{

template<typename DDT>
class Vertex_const_iterator : public std::iterator<std::input_iterator_tag, Vertex_const_iterator<DDT> >
{
public:
    typedef typename DDT::Traits                     Traits;
    typedef typename DDT::Id                     Id;
    typedef typename DDT::Tile_vertex_const_iterator Tile_vertex_const_iterator;
    typedef typename DDT::Tile_const_iterator        Tile_const_iterator;
    typedef typename DDT::Point                      Point;
    typedef typename Traits::Flag_V                    Flag_V;
    typedef typename Traits::Data_V                    Data_V;

private:
    Tile_const_iterator begin_;
    Tile_const_iterator end_;
    Tile_const_iterator tile_;
    Tile_vertex_const_iterator vertex_;

public:
    Vertex_const_iterator(Tile_const_iterator begin, Tile_const_iterator end)
        : begin_(begin), end_(end), tile_(begin), vertex_()
    {
        if(tile_ != end_)
        {
            vertex_ = tile_->vertices_begin();
            advance_to_main();
        }
    }

    Vertex_const_iterator(Tile_const_iterator begin, Tile_const_iterator end, Tile_const_iterator tile)
        : begin_(begin), end_(end), tile_(tile), vertex_()
    {
        if(tile_ != end_)
        {
            vertex_ = tile_->vertices_begin();
            advance_to_main();
        }
    }

    Vertex_const_iterator(Tile_const_iterator begin, Tile_const_iterator end, Tile_const_iterator tile, Tile_vertex_const_iterator vertex)
        : begin_(begin), end_(end), tile_(tile), vertex_(vertex)
    {
        // do not enforce main here !
    }

    Vertex_const_iterator(const Vertex_const_iterator& v)
        : begin_(v.begin_), end_(v.end_), tile_(v.tile_), vertex_(v.vertex_)
    {
        // do not enforce main here !
    }

    Vertex_const_iterator& advance_to_main()
    {
        while(tile_ != end_)
        {
            if(vertex_ == tile_->vertices_end())
            {
                if (++tile_ != end_) vertex_ = tile_->vertices_begin();
            }
            else if(tile_->vertex_is_main(vertex_))
            {
                break;
            }
            else
            {
                ++vertex_;
            }
        }
        return *this;
    }

    Vertex_const_iterator& operator++()
    {
        assert(tile_ != end_);
        ++vertex_;
        return advance_to_main();
    }

    Vertex_const_iterator operator++(int)
    {
        Vertex_const_iterator tmp(*this);
        operator++();
        return tmp;
    }

    Vertex_const_iterator& operator+=(int n)
    {
        assert(tile_ != end_);
        for(auto vit = tile_->vertices_begin(); vit != vertex_; ++vit)
            if(tile_->vertex_is_main(vit))
                ++n;
        int num_main_vertices = tile_->number_of_main_vertices();
        while(n >= num_main_vertices)
        {
            n -= num_main_vertices;
            ++tile_;
            num_main_vertices = tile_->number_of_main_vertices();
        }
        vertex_ = tile_->vertices_begin();
        advance_to_main();
        for(; n>0 ; --n)
            ++(*this);
        return *this;
    }

    bool operator==(const Vertex_const_iterator& rhs) const
    {
        return begin_ == rhs.begin_
               && tile_ == rhs.tile_
               && end_ == rhs.end_
               && (tile_ == end_ || vertex_==rhs.vertex_);
    }

    bool operator<(const Vertex_const_iterator& rhs) const
    {
        return vertex_ < rhs.vertex();
    }

    bool operator!=(const Vertex_const_iterator& rhs) const { return !(*this == rhs); }
    Vertex_const_iterator& operator*() { return *this; }
    Vertex_const_iterator* operator->() { return this; }

    Vertex_const_iterator main() const
    {
        assert(tile_ != end_);
        if(is_main()) return *this;
        Id id = main_id();
        for(Tile_const_iterator tile = begin_; tile != end_; ++tile)
            if (tile->id() == id)
                return Vertex_const_iterator(begin_, end_, tile, tile->locate_vertex(*tile_, vertex_));
        return Vertex_const_iterator(begin_, end_, end_);
    }

    Id main_id() const
    {
        assert(tile_ != end_);
        return tile_->id(vertex_);
    }

    Id gid() const
    {
        assert(tile_ != end_);
        return tile_->gid(vertex_);
    }

    Point point() const { return vertex_->point(); }

    Tile_const_iterator        tile  () const { return tile_;   }
    Tile_vertex_const_iterator vertex() const { return vertex_; }

    bool is_local()    const { return tile_->vertex_is_local(vertex_); }
    bool is_foreign()  const { return tile_->vertex_is_foreign(vertex_); }
    bool is_main()     const { return tile_->vertex_is_main(vertex_); }
    bool is_infinite() const { return tile_->vertex_is_infinite(vertex_); }

    const  Data_V & vertex_data() const
    {
        return tile_->datav(vertex_);
    }


    bool is_valid()    const
    {
        return tile_ == end_ || vertex_ != tile_->vertices_end();
    }
};

}

#endif // DDT_VERTEX_CONST_ITERATOR_HPP
