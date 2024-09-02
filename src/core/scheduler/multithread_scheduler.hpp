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
#ifndef DDT_MULTITHREAD_SCHEDULER_HPP
#define DDT_MULTITHREAD_SCHEDULER_HPP

#include <map>
#include <vector>
#include <chrono>
#include <unordered_set>
#include "thread_pool.hpp"

namespace ddt
{

template<typename Tile_>
struct multithread_scheduler
{
    typedef Tile_ Tile;
    typedef typename Tile::Vertex_const_handle_and_id Vertex_const_handle_and_id;
    typedef typename Tile::Vertex_const_handle Vertex_const_handle;
    typedef typename Tile::Vertex_handle Vertex_handle;
    typedef typename Tile::Cell_const_handle Cell_const_handle;
    typedef typename Tile::Point_and_id Point_and_id;
    typedef typename Tile::Point Point;
    typedef typename Tile::Id Id;
    multithread_scheduler(int n_threads=0) : pool(n_threads), timeout_(60)
    {
        pool.init();
    }
    template<class Duration> void timeout(Duration t)
    {
        timeout_ = t;
    }
    inline int number_of_threads() const
    {
        return pool.number_of_threads();
    }
    ~multithread_scheduler()
    {
        pool.shutdown();
    }

    void send(const Point& p, Id id, Id target)
    {
        inbox[target].emplace_back(p,id);
    }

    void send(const Point& p, Id id)
    {
        send(p,id,id);
    }

    std::function<int(Tile&)>
    insert(bool do_simplify)
    {
        return [this, do_simplify](Tile& tile)
        {
            std::vector<Point_and_id> received;
            inbox[tile.id()].swap(received);
            return int(tile.insert(received, do_simplify));
        };
    }

    template<typename F>
    std::function<int(Tile&)>
    splay_func(F&& f, bool do_simplify)
    {
        return [this,f,do_simplify](Tile& tile)
        {
            std::vector<Point_and_id> received;
            inbox[tile.id()].swap(received);
            if(!tile.insert(received, do_simplify)) return 0;
            std::vector<Vertex_const_handle_and_id> vertices;
            f(tile, std::back_inserter(vertices));
            std::map<Id, std::vector<Point_and_id>> outgoing;
            int count = tile.send_one(outgoing, vertices);
            for(auto& p : outgoing) inbox[p.first].append(p.second);
            return count;
        };
    }

    template<typename F>
    std::function<int(Tile&)>
    insert_simplified(F&& f, bool do_simplify)
    {
        return [this,f,do_simplify](Tile& tile)
        {
            std::vector<Point_and_id> received;
            inbox[tile.id()].swap(received);
            std::vector<Vertex_handle> inserted;
            tile.insert_simplified(received, inserted);
            if(inserted.empty()) return 0;
            //int inserted_size = inserted.size();
            std::unordered_set<Vertex_handle> removed;
            if(do_simplify)
            {
                tile.simplify_range(inserted.begin(), inserted.end(), removed);
                inserted.erase(std::remove_if(inserted.begin(), inserted.end(), [&](Vertex_handle v)
                {
                    return removed.find(v) != removed.end();
                } ), inserted.end());
            }
            std::map<Id, std::unordered_set<Vertex_const_handle>> vertices;
            f(tile, inserted.begin(), inserted.end(), vertices);
            int count = 0;
            for(auto& v : vertices)
            {
                std::vector<Point_and_id> outgoing;
                count += tile.send_one(v.first, v.second.begin(), v.second.end(), std::back_inserter(outgoing));
                inbox[v.first].append(outgoing);
            }
            return count;
        };
    }

    std::function<int(Tile&)>
    insert_splay(bool do_simplify)
    {
        return [this,do_simplify](Tile& tile)
        {
            std::vector<Point_and_id> received;
            inbox[tile.id()].swap(received);
            std::map<Id, std::vector<Vertex_const_handle>> vertices;
            if(!tile.insert_splay(received, vertices, do_simplify)) return 0;
            int count = 0;
            for(auto& v : vertices)
            {
                std::vector<Point_and_id> outgoing;
                count += tile.send_one(v.first, v.second.begin(), v.second.end(), std::back_inserter(outgoing));
                inbox[v.first].append(outgoing);
            }
            return count;
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
            std::map<Id, std::vector<Point_and_id>> outgoing;
            int count = tile.send_all(outgoing, vertices, begin, end);
            for(auto& p : outgoing) inbox[p.first].append(p.second);
            return count;
        };
    }

    template<class InputIt1, class InputIt2, class BinaryOp>
    void for_each(InputIt1 first1, InputIt1 last1, InputIt2 first2, BinaryOp binop)
    {
        std::vector<std::future<void>> futures;
        while(first1 != last1)
            futures.push_back(pool.submit(binop, std::ref(*first1++), std::ref(*first2++)));
        for(auto& f: futures) f.wait();
    }

    template<class InputIt, class UnaryFunction>
    void for_each(InputIt first, InputIt last, UnaryFunction f)
    {
        std::vector<std::future<void>> futures;
        while(first != last)
            futures.push_back(pool.submit(f, std::ref(*first++)));
        for(auto& f: futures) f.wait();
    }

    template<class InputIt, class UnaryPredicate>
    bool all_of(InputIt first, InputIt last, UnaryPredicate p)
    {
        std::vector<std::future<bool>> futures;
        while(first != last)
            futures.push_back(pool.submit(p, std::ref(*first++)));
        bool result = true;
        for(auto& f: futures) result = result && f.get();
        return result;
    }

    template<class InputIt, class T, class UnaryOp>
    T transform_sum(InputIt first, InputIt last,
                    T init, UnaryOp unary_op)
    {
        return transform_reduce(first, last, init, std::plus<T>(), unary_op);
    }

    template<class InputIt, class UnaryOp>
    typename std::result_of<UnaryOp(typename InputIt::reference)>::type
    transform_sum(InputIt first, InputIt last, UnaryOp unary_op)
    {
        typedef typename std::result_of<UnaryOp(typename InputIt::reference)>::type T;
        T init {};
        return transform_reduce(first, last, init, std::plus<T>(), unary_op);
    }

    template<class InputIt, class T, class BinaryOp, class UnaryOp>
    T transform_reduce(InputIt first, InputIt last,
                       T init, BinaryOp binop, UnaryOp unary_op)
    {
        std::vector<std::future<T>> futures;
        while(first != last)
            futures.push_back(pool.submit(unary_op, std::ref(*first++)));
        for(auto& f: futures) init = binop(init, f.get());
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
        result_type count = 0;
        if (begin == end) return count;
        int itcount, checkcount;
        do
        {
            std::map<Id, std::future<int>> futures;
            for(InputIt it = begin; it != end; ++it)
            {
                futures[it->id()] = pool.submit(unary_op, std::ref(*it));
            }
            bool loop = true;
            InputIt it = begin, itend = begin;
            while(loop || it != itend)
            {
                loop = false;
                if (futures[it->id()].wait_for(timeout_) != std::future_status::ready)
                {
                    itend = it;
                    if (++itend == end) itend = begin;
                    loop = true;
                }
                else
                {
                    itcount = futures[it->id()].get();
                    if (itcount)
                    {
                        count += itcount;
                        itend = it;
                        if (++itend == end) itend = begin;
                        loop = true;
                    }
                    futures[it->id()] = pool.submit(unary_op, std::ref(*it));
                }
                if (++it == end) it = begin;
            }
            checkcount = 0;
            for(InputIt it = begin; it != end; ++it)
            {
                checkcount += futures[it->id()].get();
            }
            count += checkcount;
        }
        while (checkcount);
        return count;
    }
private:
    std::map<Id, safe<std::vector<Point_and_id>>> inbox;
    thread_pool pool;
    std::chrono::milliseconds timeout_;
};

}

#endif // DDT_MULTITHREAD_SCHEDULER_HPP
