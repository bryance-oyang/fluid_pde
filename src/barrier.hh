#ifndef BARRIER_H
#define BARRIER_H

#include <condition_variable>
#include <mutex>

class ThreadBarrier {
public:
	std::mutex mutex;
	std::condition_variable cond;
	int gate_id;
	int nwaiting;
	const int nthread;

	ThreadBarrier(int nthread) : gate_id{0}, nwaiting{0}, nthread{nthread} {};

	void wait() {
		std::unique_lock<std::mutex> lock{mutex};

		nwaiting++;
		if (nwaiting == nthread) {
			nwaiting = 0;
			gate_id++;
			cond.notify_all();
		} else {
			int current_gate = gate_id;
			do {
				cond.wait(lock);
			} while (current_gate == gate_id);
		}
	}
};

#endif /* BARRIER_H */
