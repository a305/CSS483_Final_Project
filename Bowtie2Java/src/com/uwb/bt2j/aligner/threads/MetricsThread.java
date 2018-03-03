package com.uwb.bt2j.aligner.threads;

import com.uwb.bt2j.aligner.ifndef;

public class MetricsThread implements Runnable {

	@Override
	public void run() {
		Timer _t(cerr, "Multiseed full-index search: ", timing);

#ifndef _WIN32
		int pid = 0;
		if(thread_stealing) {
			pid = getpid();
			write_pid(thread_stealing_dir.c_str(), pid);
			thread_counter = 0;
		}
#endif
		
		for(int i = 0; i < nthreads; i++) {
			tids.push_back(i);
#ifdef WITH_TBB
			tps[i].tid = i;
			tps[i].done = &all_threads_done;

			if(bowtie2p5) {
				threads.push_back(new std::thread(multiseedSearchWorker_2p5, (void*)&tps[i]));
			} else {
				threads.push_back(new std::thread(multiseedSearchWorker, (void*)&tps[i]));
			}
			threads[i]->detach();
			SLEEP(10);
#else
			// Thread IDs start at 1
			if(bowtie2p5) {
				threads.push_back(new tthread::thread(multiseedSearchWorker_2p5, (void*)&tids.back()));
			} else {
				threads.push_back(new tthread::thread(multiseedSearchWorker, (void*)&tids.back()));
			}
#endif
		}

#ifndef _WIN32
		if(thread_stealing) {
			int orig_threads = nthreads;
			thread_monitor(pid, orig_threads, tids, threads);
		}
#endif
	
#ifdef WITH_TBB
		while(all_threads_done < nthreads) {
			SLEEP(10);
		}
#else
		for (int i = 0; i < nthreads; i++) {
			threads[i]->join();
		}
#endif
		for (int i = 0; i < nthreads; ++i) {
			delete threads[i];
		}

#ifndef _WIN32
		if(thread_stealing) {
			del_pid(thread_stealing_dir.c_str(), pid);
		}
#endif
	}
}
