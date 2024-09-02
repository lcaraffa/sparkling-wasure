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
#ifndef DDT_SEQUENTIAL_SCHEDULER_HPP
#define DDT_SEQUENTIAL_SCHEDULER_HPP

#include <map>
#include <vector>

namespace ddt
{

template<typename Tile_>
struct sequential_scheduler
{
    typedef Tile_ Tile;
    typedef typename Tile::Vertex_const_handle_and_id Vertex_const_handle_and_id;
    typedef typename Tile::Vertex_const_handle Vertex_const_handle;
    typedef typename Tile::Point_id Point_id;
    typedef typename Tile::Point Point;
    typedef typename Tile::Id Id;

    sequential_scheduler(int /*unused*/ = 0) {}
    inline int number_of_threads() const
    {
        return 0;
    }

    std::function<int(Tile&)>
    insert(bool do_simplify)
    {
        return [this, do_simplify](Tile& tile)
        {
            std::vector<Point_id> received;
            inbox[tile.id()].swap(received);
            return int(tile.insert(received, do_simplify));
        };
    }

    std::function<int(Tile&)>
    insert_splay(bool do_simplify)
    {
        std::cerr << "sequential_scheduler::insert_splay defaults to sequential_scheduler::insert" << std::endl;
        return insert(do_simplify);
    }

    template<typename F>
    std::function<int(Tile&)>
    insert_simplified(F&& f, bool do_simplify)
    {
        std::cerr << "sequential_scheduler::insert_simplified defaults to sequential_scheduler::insert" << std::endl;
        return insert(do_simplify);
    }

    template<typename F>
    std::function<int(Tile&)>
    splay_func(F&& f, bool do_simplify)
    {
        return [this,f,do_simplify](Tile& tile)
        {
            std::vector<Point_id> received;
            inbox[tile.id()].swap(received);
            if(!tile.insert(received, do_simplify)) return 0;
            std::vector<Vertex_const_handle_and_id> outgoing;
            f(tile, std::back_inserter(outgoing));
            return tile.send_one(inbox, outgoing);
        };
    }

    template<typename Id_iterator, typename F, typename... Args>
    std::function<int(const Tile&)>
    send_all_func(Id_iterator begin, Id_iterator end, F&& f, Args&&... args)
    {
        return [this,f,begin,end,args...](const Tile& tile)
        {
            std::vector<Vertex_const_handle> vertices;
            f(tile, std::back_inserter(vertices), args...);
            return tile.send_all(inbox, vertices, begin, end);
        };
    }

    template<class InputIt, class UnaryFunction>
    void for_each(InputIt first, InputIt last, UnaryFunction f) const
    {
        std::for_each(first, last, f);
    }

    template<class InputIt, class UnaryPredicate>
    bool all_of(InputIt first, InputIt last, UnaryPredicate p) const
    {
        return std::all_of(first, last, p);
    }

    template<class InputIt, class T, class UnaryOp>
    T transform_sum(InputIt first, InputIt last,
                    T init, UnaryOp unary_op) const
    {
        return transform_reduce(first, last, init, std::plus<T>(), unary_op);
    }

    template<class InputIt, class UnaryOp>
    typename std::result_of<UnaryOp(typename InputIt::reference)>::type transform_sum(InputIt first, InputIt last, UnaryOp unary_op) const
    {
        typedef typename std::result_of<UnaryOp(typename InputIt::reference)>::type T;
        T init {};
        return transform_reduce(first, last, init, std::plus<T>(), unary_op);
    }

    template<class InputIt, class T, class BinaryOp, class UnaryOp>
    T transform_reduce(InputIt first, InputIt last,
                       T init, BinaryOp binop, UnaryOp unary_op) const
    {
        while (first != last) init = binop(init, unary_op(*first++));
        return init;
    }

    template<typename InputIt, typename UnaryOp>
    typename std::result_of<UnaryOp(typename InputIt::reference)>::type
    for_each_rec(InputIt begin, InputIt end, UnaryOp unary_op)
    {
        typedef typename std::result_of<UnaryOp(typename InputIt::reference)>::type result_type;
        result_type count = transform_sum(begin, end, unary_op), c;
        while((c = transform_sum(begin, end, unary_op)))
            count += c;
        return count;
    }

    template<typename InputIt, typename UnaryOp>
    typename std::result_of<UnaryOp(typename InputIt::reference)>::type
    for_each_rec2(InputIt begin, InputIt end, UnaryOp unary_op)
    {
        typedef typename std::result_of<UnaryOp(typename InputIt::reference)>::type result_type;
        result_type count = transform_sum(begin, end, unary_op), c;
        InputIt itend = end;
        for(InputIt it = begin; it != itend; ++it)
        {
            if (it == end) it = begin;
            if((c = unary_op(*it)))
            {
                count += c;
                itend = it;
            }
        }
        return count;
    }

    void send(const Point& p, Id id, Id target)
    {
        inbox[target].emplace_back(p,id);
    }

    void send(const Point& p, Id id)
    {
        send(p,id,id);
    }
private:
    std::map<Id, std::vector<Point_id>> inbox;
};

}

#endif // DDT_SEQUENTIAL_SCHEDULER_HPP
