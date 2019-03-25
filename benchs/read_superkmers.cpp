
#include <cstddef>
#include <iostream>
#include <cstring>
#include <thread>
#include <algorithm>

#include "gatbl/sys/file.hpp"
#include "gatbl/fastx.hpp"


using namespace gatbl;

template<typename F>
inline void
process_fastq(std::string filename, F&& callback)
{
    sys::file_descriptor fd(filename);
    auto data = fd.mmap<const char>();
    // FIXME: those are optional but we really want to knwow if they work for testing
    assert(data.advise_hugepage());
    assert(data.advise_sequential());
    for (auto& rec : seq_record_subrange<fastq_record<const char*>>(data)) {
        callback(rec);
    }
}

using namespace std;

typedef uint64_t kmer_type;

inline static char char_to_nt(char c)
{
    return (c>>1) & 3;
}

inline static char char_to_nt_rc(char c)
{
    return ((c>>1)^2) & 3;
}

// to avoid being optimized out
template<typename T>
void
use(T&& t)
{
        __asm__ __volatile__("" ::"g"(t));
}


void emit_kmer(kmer_type kmer)
{
    // do what we want with that kmer
    //std::cout << kmer << std::endl;
    use(kmer);
}

bool is_minimizer_allowed(uint32_t mmer)
{
    // TODO
    // for now is strict lexical order, not KMC's order
    return true;
}

inline static uint64_t revcomp (const uint64_t & x, int k)
{
    uint64_t res = x;
    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2*(32-k))) ;
}

void init_minimizers(int _minimizerSize, vector<uint32_t> &_mmer_lut)
{
    uint32_t _mask  = (1 << (2*_minimizerSize)) - 1;

    uint32_t nbminims_total = (1 << (2*_minimizerSize));

    for(uint32_t ii=0; ii< nbminims_total; ii++)
    {
        uint32_t mmer = ii;
        uint32_t rev_mmer = revcomp(mmer, _minimizerSize);
        if(rev_mmer < mmer) mmer = rev_mmer;

        if (!is_minimizer_allowed(mmer),_minimizerSize)
            mmer = _mask;

        _mmer_lut[ii] = mmer;
    }
}

void compute_minimizer(uint64_t kmer, uint32_t &minimizer, unsigned int& minimizer_position, unsigned int k, const vector<uint32_t> &_mmer_lut, unsigned int _minimizerSize)
{
    minimizer = 0;
    uint32_t mmer;
    uint32_t  mmerMask  = (1 << (2*_minimizerSize)) - 1;
    for (size_t i=0; i<k-_minimizerSize+1; ++i)
    {
        mmer    = _mmer_lut[(kmer >>i)& mmerMask];
        if (mmer < minimizer)
        {
            minimizer = mmer;
            minimizer_position = i;
        }
    }
}

void read_to_superkmer(const char* read_start, const char* read_end,
        const vector<uint32_t> &_mmer_lut, unsigned int k, unsigned int _minimizerSize, 
        uint64_t &nb_kmers_thread, uint64_t &nb_superkmers_thread)
{
    kmer_type kmer = 0, kmer_rc = 0;
    kmer_type kmerMask = (1LL << (k*2)) - 1;
    uint32_t  mmerMask  = (1 << (2*_minimizerSize)) - 1;
    unsigned int minimizer_position;
    uint32_t minimizer, mmer;
    string superkmer = "";
    uint64_t last_kmer = 0;

    size_t buffer_index = 0;
    size_t readlen = read_end - read_start;
    while (buffer_index < readlen-k+1)
    {
        bool found_N = false;
        // prime the first kmer
        for (size_t i=buffer_index; i<buffer_index+k; ++i)
        {
            if (unlikely(read_start[i] == 'N')) {found_N = true; buffer_index = i+1; break;}
            unsigned char c  = char_to_nt(read_start[i]);
            kmer    = (kmer<<2) + c;
            kmer_rc = (kmer>>2) + (c^2);
        }
        if (unlikely(found_N)) continue;

        emit_kmer(std::min(kmer,kmer_rc));
        nb_kmers_thread++;
        compute_minimizer(kmer, minimizer, minimizer_position, k, _mmer_lut, _minimizerSize);

        buffer_index++;

        // iterate next kmers of that read
        while (buffer_index < readlen-k+1)
        {
            size_t i = buffer_index+k;
            if (unlikely(read_start[i] == 'N')) { buffer_index = i+1; break;} // no need to set found_N here, we'll continue 
            unsigned char c  = char_to_nt(read_start[i]);
            kmer    = (( kmer << 2) +  c) & kmerMask;
            kmer_rc = (( kmer >> 2) +  (c^2)) & kmerMask;
            superkmer += c;
            emit_kmer(std::min(kmer,kmer_rc));
            nb_kmers_thread++;

            mmer    = _mmer_lut[kmer & mmerMask];
            minimizer_position--;
            if (mmer <= minimizer)
            {
                if (mmer != minimizer) 
                {
                    use(superkmer); // emit the superkmer here
                    nb_superkmers_thread++;
                    minimizer = mmer;
                    superkmer = last_kmer; // it works?! TODO check
                }
                minimizer_position = i-_minimizerSize+1;
            }
            else
            {
                if (minimizer_position < 0)
                {
                    uint32_t old_minimizer = minimizer;
                    compute_minimizer(kmer, minimizer, minimizer_position, k, _mmer_lut, _minimizerSize);
                    if (minimizer != old_minimizer)
                    {
                        use(superkmer); // emit the superkmer here
                        nb_superkmers_thread++;
                        superkmer = last_kmer; // it works?! TODO check
                    }
                }
            }
            last_kmer = ((last_kmer<<2) + c ) & kmerMask;
            buffer_index++;
        }

    }

    use(superkmer); // emit the last superkmer
    nb_superkmers_thread++;
}

int
main(int argc, char** argv)
{
    if (argc < 2)
        return 1;

    int nb_threads = 6;

    if (argc == 3)
        nb_threads = atoi(argv[2]);

    int _minimizerSize = 8;
    int k = 25;


    vector<vector<const char*>> read_positions(nb_threads); // fill a [for now non-circular] buffer of read position, written by the IO thread and read by the worker threads
    vector<unsigned int> nb_kmers_threads(nb_threads), nb_superkmers_threads(nb_threads);
    uint32_t nbminims_total = (1 << (2*_minimizerSize));
    vector<uint32_t> _mmer_lut(nbminims_total); 
    init_minimizers(_minimizerSize, _mmer_lut);
    uint64_t nb_reads = 0;
    vector<uint64_t> nb_read_positions(nb_threads);
    vector<std::thread> threads;
    vector<bool> thread_over(nb_threads);
    vector<bool> thread_caught_up(nb_threads);
    const uint64_t read_positions_buffer_size = (1<<26); // process 70M reads/thread at a time (memory usage of 140MB/thread to remember parsed read positions)

    for (int thread_id = 0; thread_id < nb_threads; thread_id++)
    { 
        thread_over[thread_id] = false;
        thread_caught_up[thread_id] = true;
        read_positions[thread_id].reserve(read_positions_buffer_size);
        nb_read_positions[thread_id] = 0;

        auto thread_func = [&thread_over, &thread_caught_up, read_positions_buffer_size, nb_threads, thread_id, &read_positions, &nb_read_positions, &nb_kmers_threads, &nb_superkmers_threads, &_mmer_lut, &k, &_minimizerSize]()
        {
            uint64_t position_index = 0, nb_kmers_thread = 0, nb_superkmers_thread = 0;
            while ((!thread_over[thread_id]) || (position_index < nb_read_positions[thread_id]))
            {
                while ((position_index >= nb_read_positions[thread_id]) && (!thread_over[thread_id])) { 
                    if (thread_caught_up[thread_id] == false) // signal that we're going to reset read_positions and nb_read_positions
                    {
                        if (!(position_index >= nb_read_positions[thread_id])) // needed to avoid false alerts (which result in slightly 
                                                                               // less sequences processed: to reproduce, uncomment 
                                                                               // that "if" statement and set read_positions_buffer_size 
                                                                               // to 1<<11), and make sure to sleep(0.1) just before the 
                                                                               // parent "if" (not after)
                            break;
                        read_positions[thread_id].clear();
                        nb_read_positions[thread_id] = 0;
                        position_index = 0;
                        thread_caught_up[thread_id] = true;
                    }
                    sleep(0.1);   // FIXME: code loops forever if that sleep() is missing; would love to figure out why; to reproduce: 
                                  // comment that line and do benchs/read_superkmers benchs/frag_1.fastq.20klines  3
                    continue;
                } // waiting on IO, shouldn't happen so often
                if (position_index >= nb_read_positions[thread_id]) // avoid case when thread is over but position_index is already too high
                    continue;
                const char* read_start = read_positions[thread_id][position_index];
                const char* read_end = read_positions[thread_id][position_index+1];
                position_index += 2;
                read_to_superkmer(read_start, read_end, _mmer_lut, k, _minimizerSize,  nb_kmers_thread, nb_superkmers_thread);
            }

            // eliminates false sharing instead of incrementing nb_kmers_threads online
            nb_kmers_threads[thread_id] = nb_kmers_thread;
            nb_superkmers_threads[thread_id] = nb_superkmers_thread;
        };

        // spawn thread
        threads.push_back(std::thread(thread_func));
    }

    process_fastq(argv[1], [nb_threads, &nb_reads, &nb_read_positions, &read_positions, &thread_caught_up](fastq_record<>& rec) { 

            int thread_id = nb_reads % nb_threads;
            read_positions[thread_id].push_back(rec.sequence().begin());
            read_positions[thread_id].push_back(rec.sequence().end());
            nb_reads ++;
            nb_read_positions[thread_id] += 2;

            // when one of the threads' buffer of read positions is about to overflow
            if (unlikely((nb_read_positions[thread_id] > 0) && (nb_read_positions[thread_id] % read_positions_buffer_size == 0)))
            {
                // signal to the threads that they need to catch up
                std::fill(thread_caught_up.begin(), thread_caught_up.end(), false);
                // wait for all threads to have caught up
                while (! std::all_of(thread_caught_up.begin(), thread_caught_up.end(), [](bool v) { return v; }))
                {
                    continue;
                }
            }
    });
    std::cout << "done parsing fastq" << std::endl;
   
    for (int i = 0; i < nb_threads; i++)
        thread_over[i] = true;
    
    uint64_t nb_kmers = 0, nb_superkmers = 0;
    for (int i = 0; i < nb_threads; i++)
    {
        threads[i].join();
        nb_kmers += nb_kmers_threads[i] ; nb_superkmers += nb_superkmers_threads[i];
    }

    std::cout << "nb kmers:      " << nb_kmers      << std::endl;
    std::cout << "nb superkmers: " << nb_superkmers << std::endl;
}
