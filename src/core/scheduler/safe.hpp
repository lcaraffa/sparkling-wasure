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
#ifndef SAFE_HP
#define SAFE_HPP

#include <mutex>

// Thread safe implementation of a Queue using a std::queue
template <typename Container>
class safe
{
private:
    Container m_container;
    std::mutex m_mutex;
public:
    typedef typename Container::value_type value_type;
    typedef typename Container::size_type size_type;
    safe() {}

    safe(const safe& other) = delete;

    ~safe() {}

    bool empty()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return m_container.empty();
    }

    size_type size()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return m_container.size();
    }

    void enqueue(value_type& t)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_container.push(t);
    }

    bool dequeue(value_type& t)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        if (m_container.empty()) return false;
        t = std::move(m_container.front());
        m_container.pop();
        return true;
    }

    template< class... Args >
    void emplace_back( Args&&... args )
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_container.emplace_back(std::forward<Args>(args)...);
    }

    void swap( Container& container )
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_container.swap(container);
    }

    void append( const Container& container )
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_container.insert(m_container.end(), container.begin(), container.end());
    }
};

#endif // SAFE_HPP
