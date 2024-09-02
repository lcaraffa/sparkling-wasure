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
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>
#include <queue>

#include "safe.hpp"

class thread_pool
{
private:
    class thread_worker
    {
    private:
        thread_pool * m_pool;
        int m_id;
    public:
        thread_worker(thread_pool * pool, const int id)
            : m_pool(pool), m_id(id)
        {
        }

        void operator()()
        {
            std::function<void()> func;
            bool dequeued;
            while (!m_pool->m_shutdown)
            {
                {
                    std::unique_lock<std::mutex> lock(m_pool->m_conditional_mutex);
                    if (m_pool->m_queue.empty())
                    {
                        m_pool->m_conditional_lock.wait(lock);
                    }
                    dequeued = m_pool->m_queue.dequeue(func);
                }
                if (dequeued)
                {
                    func();
                }
            }
        }
    };

    std::vector<std::thread> m_threads;
    bool m_shutdown;
    safe<std::queue<std::function<void()>>> m_queue;
    std::mutex m_conditional_mutex;
    std::condition_variable m_conditional_lock;
public:
    thread_pool(const int n_threads)
        : m_threads(n_threads?n_threads:std::thread::hardware_concurrency()), m_shutdown(false)
    {
    }

    thread_pool(const thread_pool &) = delete;
    thread_pool(thread_pool &&) = delete;

    thread_pool & operator=(const thread_pool &) = delete;
    thread_pool & operator=(thread_pool &&) = delete;

    inline int number_of_threads() const
    {
        return m_threads.size();
    }

    // Inits thread pool
    void init()
    {
        for (size_t i = 0; i < m_threads.size(); ++i)
        {
            m_threads[i] = std::thread(thread_worker(this, i));
        }
    }

    // Waits until threads finish their current task and shutdowns the pool
    void shutdown()
    {
        m_shutdown = true;
        m_conditional_lock.notify_all();
        for (size_t i = 0; i < m_threads.size(); ++i)
        {
            m_threads[i].join();
        }
    }

    // Submit a function to be executed asynchronously by the pool
    template<typename F, typename...Args>
    auto submit(F&& f, Args&&... args) -> std::future<decltype(f(args...))>
    {
        // Create a function with bounded parameters ready to execute
        std::function<decltype(f(args...))()> func = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
        // Encapsulate it into a shared ptr in order to be able to copy construct / assign
        auto task_ptr = std::make_shared<std::packaged_task<decltype(f(args...))()>>(func);
        // Wrap packaged task into void function
        std::function<void()> wrapper_func = [task_ptr]()
        {
            (*task_ptr)();
        };
        // Enqueue generic wrapper function
        m_queue.enqueue(wrapper_func);
        // Wake up one thread if its waiting
        m_conditional_lock.notify_one();
        // Return future from promise
        return task_ptr->get_future();
    }
};

#endif // THREAD_POOL_HPP
