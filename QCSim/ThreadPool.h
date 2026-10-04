#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "QubitRegisterCalculator.h"

namespace QC
{
template <typename Return = double> class ThreadPool
{
  public:
    explicit ThreadPool(size_t numThreads = 0)
    {
        if (numThreads == 0)
            numThreads = QubitRegisterCalculator<>::GetNumberOfThreads();
        if (numThreads == 0)
            numThreads = 4;

        try
        {
            for (size_t i = 0; i < numThreads; ++i)
            {
                workers.emplace_back([this] {
                    currentPool = this;
                    for (;;)
                    {
                        std::function<void()> task;
                        {
                            std::unique_lock<std::mutex> lock(mtx);
                            cv.wait(lock, [this] { return stopped || !tasks.empty(); });
                            if (stopped && tasks.empty())
                                return;
                            // taskRunning = true;
                            task = std::move(tasks.front());
                            tasks.pop();
                        }
                        task();
                        /*
                        {
                            std::unique_lock<std::mutex> lock(mtx);
                            taskRunning = false;
                        }
                        cvDone.notify_all();
                        */
                    }
                });
            }
        }
        catch (...)
        {
            {
                std::lock_guard<std::mutex> lock(mtx);
                stopped = true;
            }
            cv.notify_all();
            for (auto &worker : workers)
                worker.join();
            throw;
        }
    }

    ~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(mtx);
            stopped = true;
        }
        cv.notify_all();
        for (auto &worker : workers)
            worker.join();
    }

    ThreadPool(const ThreadPool &) = delete;
    ThreadPool &operator=(const ThreadPool &) = delete;
    ThreadPool(ThreadPool &&) = delete;
    ThreadPool &operator=(ThreadPool &&) = delete;

    template <typename F> std::future<Return> Enqueue(F &&f)
    {
        auto task = std::make_shared<std::packaged_task<Return()>>(std::forward<F>(f));
        std::future<Return> theFuture = task->get_future();
        {
            std::unique_lock<std::mutex> lock(mtx);
            tasks.emplace([task = std::move(task)]() { (*task)(); });
        }
        cv.notify_one();
        return theFuture;
    }

    /*
    bool HasWork()
    {
        std::unique_lock<std::mutex> lock(mtx);
        return !tasks.empty() || taskRunning;
    }

    bool WaitForFinish()
    {
        std::unique_lock<std::mutex> lock(mtx);
        if (tasks.empty() && !taskRunning)
            return true;
        cvDone.wait(lock, [this] { return tasks.empty() && !taskRunning; });
        return tasks.empty() && !taskRunning;
    }
    */

    size_t GetThreadCount() const
    {
        return workers.size();
    }
    bool IsWorkerThread() const
    {
        return currentPool == this;
    }

  private:
    inline static thread_local const ThreadPool *currentPool = nullptr;
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex mtx;
    std::condition_variable cv;
    // std::condition_variable cvDone;
    bool stopped = false;
    // bool taskRunning = false;
};
} // namespace QC
