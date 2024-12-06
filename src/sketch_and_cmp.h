#pragma once
#include "dashing.h"

// Import types from sketch and bns namespaces
using ::sketch::hll::EstimationMethod;
using ::sketch::hll::JointEstimationMethod;
using bns::EmissionFormat;
using bns::EmissionType;
using bns::Spacer;
using bns::KSeqBufferHolder;
using namespace sketch;

namespace bns {

/**
 * @brief Emits distance matrix entries to file in specified format
 * 
 * @param pairfi File descriptor to write to
 * @param ptr Pointer to distance values
 * @param hs Number of entries
 * @param index Current row index
 * @param str String buffer for output
 * @param inpaths Vector of input paths
 * @param emit_fmt Output format (TSV or upper triangular)
 * @param use_scientific Whether to use scientific notation
 * @param buffer_flush_size Size at which to flush buffer
 * @return size_t Index of current row
 */
template<typename FType=float, typename=typename std::enable_if<std::is_floating_point<FType>::value>::type>
size_t submit_emit_dists(int pairfi, const FType *ptr, u64 hs, size_t index, ks::string &str, const std::vector<std::string> &inpaths, EmissionFormat emit_fmt, bool use_scientific, const size_t buffer_flush_size=BUFFER_FLUSH_SIZE) {
    auto &strref = inpaths[index];
    str += strref;
    constexpr const char *fmt = "\t%.6g";
    if(emit_fmt == UT_TSV) {
        {
            u64 k;
            for(k = 0; k < index + 1;  ++k, kputsn_("\t-", 2, reinterpret_cast<kstring_t *>(&str)));
            for(k = 0; k < hs - index - 1; str.sprintf(fmt, ptr[k++]));
        }
    } else { // emit_fmt == UPPER_TRIANGULAR
        if(strref.size() < 9)
            str.append(9 - strref.size(), ' ');
        for(u64 k = 0; k < hs - index - 1; str.sprintf(fmt, ptr[k++]));
    }
    str.putc_('\n');
    str.flush(pairfi);
    return index;
}

// Forward declaration of distance loop function
template<typename SketchType>
void dist_loop(std::FILE *&ofp, std::string ofpname, SketchType *sketches, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int, const size_t buffer_flush_size, size_t nq);

using namespace sketch;
using namespace hll;

/**
 * @brief Converts byte size log2 to sketch size parameter
 *
 * @param nblog2 Log2 of number of bytes
 * @param sketch Type of sketch
 * @return size_t Sketch size parameter
 */
static size_t bytesl2_to_arg(int nblog2, Sketch sketch) {
    switch(sketch) {
        case HLL: case WIDE_HLL: return nblog2;
        case BLOOM_FILTER: return nblog2 + 3; // 8 bits per byte
        case RANGE_MINHASH: return size_t(1) << (nblog2 - 3); // 8 bytes per minimizer
        case COUNTING_RANGE_MINHASH: return (size_t(1) << (nblog2)) / double(sizeof(uint64_t) + sizeof(uint32_t));
        case BB_MINHASH:
            return nblog2 - std::floor(std::log2(gargs.bbnbits / 8));
        case BB_SUPERMINHASH:
            return size_t(1) << (nblog2 - int(std::log2(gargs.bbnbits / 8)));
        case FULL_KHASH_SET: return 16; // Reserve hash set size a bit. Mostly meaningless, resizing as necessary.
        case HYPERMINHASH: {
            switch(gargs.bbnbits) {
                case 8: return nblog2;
                case 16: return nblog2 - 1;
                case 32: return nblog2 - 2;
                case 64: return nblog2 - 3;
                default: {
                    if(gargs.bbnbits < 8) {gargs.bbnbits = 8; return nblog2;}
                    if(gargs.bbnbits < 16) {gargs.bbnbits = 16; return nblog2 - 1;}
                    if(gargs.bbnbits < 32) {gargs.bbnbits = 32; return nblog2 - 2;}
                    gargs.bbnbits = 64;
                    return nblog2 - 3;
                }
            }
        }
        default: {
            // Buffer for error message
            char buf[128];
            // Format error message - either "Not such sketch" if sketch type is invalid,
            // or the name of the unsupported sketch type
            std::sprintf(buf, "Sketch %s not yet supported.\n", 
                (size_t(sketch) >= (sizeof(sketch_names) / sizeof(char *)) ? 
                 "Not such sketch": sketch_names[sketch]));
            // Throw unrecoverable error with the message
            UNRECOVERABLE_ERROR(buf);
            // Return invalid value (should never be reached due to error)
            return -1337;
        }
    }

}

/**
 * @brief Performs distance calculations between sketches loaded from a binary file
 *
 * @param labels Vector of sequence/sample labels
 * @param datapath Path to binary file containing sketches
 * @param pairofp File pointer for writing pairwise distances
 * @param outpath Output path for distance matrix
 * @param k K-mer size used for sketching
 * @param estim Method for cardinality estimation
 * @param jestim Method for joint cardinality estimation
 * @param result_type Type of distance/similarity metric to compute
 * @param emit_fmt Format for emitting results (e.g. TSV, binary)
 * @param nthreads Number of threads to use
 * @param otherpath Optional path to query sketches for asymmetric comparisons
 */
template<typename SketchType>
void dist_by_seq(std::vector<std::string> &labels, std::string datapath,
                 std::FILE *pairofp, std::string outpath, int k,
                 EstimationMethod estim, JointEstimationMethod jestim, EmissionType result_type, EmissionFormat emit_fmt,
                 unsigned nthreads, std::string otherpath)
{
    // Open binary sketch file
    gzFile sfp = gzopen(datapath.data(), "rb");
    if(!sfp) throw sketch::ZlibError(std::string("Failed to open file at ") + datapath);

    // Initialize vectors to store sketches and query names
    std::vector<SketchType> sketches;
    std::vector<std::string> qnames;
    sketches.reserve(labels.size());

    // Read sketches for reference sequences
    while(sketches.size() < labels.size()) {
        sketches.emplace_back(sfp);
        set_estim_and_jestim(sketches.back(), estim, jestim);
    }

    // Handle optional query sketches
    if(otherpath.size()) {
        // Load query sequence names
        qnames = get_paths((otherpath + ".names").data());
        if(qnames.empty()) UNRECOVERABLE_ERROR("Can't compare with empty qnames");

        // Read query sketches
        sketches.reserve(qnames.size() + sketches.size());
        for(size_t i = 0; i < qnames.size(); ++i) {
            sketches.emplace_back(sfp);
            set_estim_and_jestim(sketches.back(), estim, jestim);
        }
        labels.insert(labels.end(), qnames.begin(), qnames.end());
    } else if(!is_symmetric(result_type)) {
        UNRECOVERABLE_ERROR("Can't perform asymmetric comparison without query paths");
    }

    const size_t nq = qnames.size();
    gzclose(sfp);

    // Write header information based on output format
    ks::string str;
    if(emit_fmt == UT_TSV && !nq) {
        // Write TSV header with sequence names
        size_t nq = 0;
        str.sprintf("##Names\t");
        for(size_t i = 0; i < labels.size() - nq; ++i) {
            str.sprintf("%s\t", labels[i].data());
        }
        str.back() = '\n';
        str.flush(fileno(pairofp));
    } else if(emit_fmt == UPPER_TRIANGULAR) { // emit_fmt == UPPER_TRIANGULAR
        std::fprintf(pairofp, "%zu\n", labels.size());
        std::fflush(pairofp);
    }

    // Calculate and write distances
    dist_loop<SketchType>(pairofp, outpath, sketches.data(), labels, /* use_scientific=*/ true, k, result_type, emit_fmt, nthreads, BUFFER_FLUSH_SIZE, nq);
}

/**
 * @brief Creates sketches for a set of sequences and emits size estimates
 *
 * @tparam SketchType Type of sketch to create (e.g. HyperLogLog, MinHash)
 * @param inpaths Input sequence file paths
 * @param cms Vector of counting sketches for filtering kmers
 * @param kseqs Buffer for reading sequences
 * @param ofp Output file pointer
 * @param sp Spacer for k-mer generation
 * @param ssarg Sketch size argument (log2)
 * @param mincount Minimum kmer count for inclusion
 * @param enct Encoding type (BONSAI, NTHASH, etc)
 * @param estim Estimation method for cardinality
 * @param jestim Joint estimation method
 * @param cache_sketch Whether to cache sketches to disk
 * @param emit_binary Whether to emit binary format
 * @param use_scientific Whether to use scientific notation
 * @param presketched_only Only read pre-computed sketches
 * @param nthreads Number of threads to use
 * @param suffix Suffix for sketch filenames
 * @param prefix Prefix for sketch filenames
 * @param canon Use canonical k-mers
 * @param spacing Spacing between k-mers
 */
template<typename SketchType>
void size_sketch_and_emit(std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp,
                         Spacer sp,
                         unsigned ssarg, unsigned mincount, EncodingType enct, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, bool emit_binary,
                         bool use_scientific,
                         bool presketched_only, unsigned nthreads, std::string suffix, std::string prefix, bool canon, std::string spacing)
{
    // Get final sketch type and validate it's not weighted
    using final_type = typename FinalSketch<SketchType>::final_type;
    static_assert(!sketch::wj::is_weighted_sketch<final_type>::value, "Can't have a weighted sketch be a final type");

    // Allocate buffer for sketches
    auto buf(std::make_unique<uint8_t[]>(inpaths.size() * sizeof(SketchType)));
    auto bufp = &buf[0];
    auto sketches = (SketchType *)bufp;
    const size_t npaths = inpaths.size();
    const uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);

    // Initialize sketches in parallel
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        new(sketches + i) SketchType(construct<SketchType>(sketch_size));
        set_estim_and_jestim(sketches[i], estim, jestim);
    }

    // Check if final type is same as sketch type
    static constexpr bool samesketch = std::is_same<SketchType, final_type>::value;
    final_type *final_sketches;

    // Initialize counters and parameters
    std::atomic<uint32_t> ncomplete;
    ncomplete.store(0);
    const unsigned k = sp.k_;
    const unsigned wsz = sp.w_;
    RollingHasher<uint64_t> rolling_hasher(k, canon);

    // Allocate final sketches array
    final_sketches =
        samesketch ? reinterpret_cast<final_type *>(sketches)
                   : static_cast<final_type *>(std::malloc(sizeof(*final_sketches) * npaths));

    // Process each input path in parallel
    OMP_PFOR_DYN
    for(size_t i = 0; i < npaths; ++i) {
        const std::string &path(inpaths[i]);
        auto &sketch = sketches[i];

        // Handle pre-sketched files
        if(presketched_only)  {
            CONST_IF(samesketch) {
                sketch.read(path);
                set_estim_and_jestim(sketch, estim, jestim); // HLL is the only type that needs this, and it's the same
            } else {
                new(final_sketches + i) final_type(path.data()); // Read from path
            }
        }
        // Process new sketches
        else {
            const std::string fpath(make_fname<SketchType>(path.data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix, enct));
            const bool isf = isfile(fpath);

            // Load cached sketch if available
            if(cache_sketch && isf) {
                LOG_DEBUG("Sketch found at %s with size %zu, %u\n", fpath.data(), size_t(1ull << sketch_size), sketch_size);
                CONST_IF(samesketch) {
                    sketch.read(fpath);
                    set_estim_and_jestim(sketch, estim, jestim);
                } else {
                    new(final_sketches + i) final_type(fpath);
                }
            }
            // Create new sketch
            else {
                const int tid = omp_get_thread_num();
                Encoder<score::Lex> enc(nullptr, 0, sp, nullptr, canon);

                // Process without counting filter
                if(cms.empty()) {
                    auto &h = sketch;
                    switch(enct) {
                    case BONSAI:
                        for_each_substr([&](const char *s) {enc.for_each([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, path, FNAME_SEP);
                        break;
                    case NTHASH:
                        for_each_substr([&](const char *s) {enc.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, path, FNAME_SEP);
                        break;
                    case RK: case CYCLIC: for_each_substr([&](const char *s) {rolling_hasher.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, &kseqs[tid]);}, path, FNAME_SEP);
                        break;
                    }
                }
                // Process with counting filter
                else {
                    CountingSketch &cm = cms.at(tid);
                    const auto lfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) sketch.addh(kmer);};
                    switch(enct) {
                    case BONSAI:
                        for_each_substr([&](const char *s) {enc.for_each(lfunc, s, &kseqs[tid]);}, path, FNAME_SEP);
                    break;
                    case NTHASH:
                        for_each_substr([&](const char *s) {enc.for_each_hash(lfunc, s, &kseqs[tid]);}, path, FNAME_SEP);
                    break;
                    case RK: case CYCLIC: default:
                        for_each_substr([&](const char *s) {rolling_hasher.for_each_hash(lfunc, s, &kseqs[tid]);}, path, FNAME_SEP);
                    break;
                    }
                    cm.clear();
                }

                // Handle final sketch storage
                CONST_IF(!samesketch) new(final_sketches + i) final_type(std::move(sketch));
                CONST_IF(samesketch) {
                    if(cache_sketch && !isf) sketch.write(fpath);
                } else if(cache_sketch) final_sketches[i].write(fpath);
            }
        }
        ++ncomplete; // Atomic
    } // End of parallel sketch construction loop

    // Finalize all sketches in parallel
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        try {
            sketch_finalize(final_sketches[i]); // Finalize each sketch
        } catch(const std::exception &ex) {
            std::cerr << "Failed to finalize sketch " << i << " for value " << inpaths[i] << ".\n"
            << "msg: " << ex.what() << '\n';
        }
    }

    // Free sequence buffers
    kseqs.free();

    // Calculate and store cardinality estimates
    auto fbuf = std::make_unique<float[]>(npaths);
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        // Get cardinality estimate from appropriate sketch type
        CONST_IF(samesketch) fbuf[i] = cardinality_estimate(sketches[i]);
        else                 fbuf[i] = cardinality_estimate(final_sketches[i]);
    }

    // Write cardinality estimates to output file
    if(emit_binary) {
        // Write binary format
        if(std::fwrite(&fbuf[0], sizeof(float), npaths, ofp) != npaths) {
            throw std::system_error(std::ferror(ofp), std::system_category(), "Failed to write cardinality estimates to file");
        }
    } else {
        // Write text format
        ks::string str("#Path\tSize (est.)\n");
        str.resize(BUFFER_FLUSH_SIZE);
        {
            const int fn(fileno(ofp));
            constexpr const char *scinotstr = "%s\t%0.12g\n";  // Format for scientific notation
            constexpr const char *stdnotstr = "%s\t%0.8f\n";   // Format for standard notation
            const char *const ptr = use_scientific ? scinotstr: stdnotstr;
            
            // Write each path and its cardinality estimate
            for(size_t i(0); i < npaths; ++i) {
                str.sprintf(ptr, inpaths[i].data(), fbuf[i]);
                if(str.size() >= BUFFER_FLUSH_SIZE) str.flush(fn);
            }
            str.flush(fn);
        }
    }

    // Close output file if not stdout
    if(ofp != stdout) std::fclose(ofp);

    // Clean up sketches
    CONST_IF(!samesketch) {
        OMP_PFOR
        for(size_t i = 0; i < npaths; ++i) {
            using T = typename std::decay<decltype(final_sketches[0])>::type;
            final_sketches[i].~T();  // Destruct final sketches
        }
        std::free(final_sketches);
    }

    // Clean up original sketches
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        sketches[i].~SketchType();
    }
} // size_sketch_and_emit

/**
 * @brief Performs sketching and distance comparisons between sequences
 *
 * @tparam SketchType Type of sketch to use (e.g. HLL, MinHash)
 * @param inpaths Input sequence file paths
 * @param cms Vector of counting sketches
 * @param kseqs Buffer for reading sequences
 * @param ofp Output file pointer for sketches
 * @param pairofp Output file pointer for pairwise distances
 * @param outpath Output path for distance matrix
 * @param sp Spacer for k-mer sampling
 * @param ssarg Sketch size argument (log2 of bytes)
 * @param mincount Minimum k-mer count threshold
 * @param estim Method for cardinality estimation
 * @param jestim Method for joint cardinality estimation
 * @param cache_sketch Whether to cache sketches to disk
 * @param result_type Type of distance/similarity metric
 * @param emit_fmt Format for emitting results
 * @param presketched_only Only load pre-computed sketches
 * @param nthreads Number of threads to use
 * @param use_scientific Use scientific notation in output
 * @param suffix Suffix for sketch filenames
 * @param prefix Prefix for sketch filenames
 * @param canon Use canonical k-mers
 * @param entropy_minimization Use entropy minimization
 * @param spacing K-mer sampling spacing pattern
 * @param nq Number of query sequences (default: 0)
 * @param enct Encoding type for k-mers (default: BONSAI)
 */
template<typename SketchType>
void dist_sketch_and_cmp(std::vector<std::string> &inpaths, std::vector<CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp, std::string outpath,
                         Spacer sp,
                         unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,
                         bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,
                         size_t nq=0, EncodingType enct=BONSAI)
{
    // nq -- number of queries
    //       for convenience, we will perform our comparisons (all-p) against (all-q) [remainder]
    //       and use the same guts for all portions of the process
    //       except for the final comparison and output.
    assert(nq <= inpaths.size());

    // Get final sketch type and validate it's not weighted
    using final_type = typename FinalSketch<SketchType>::final_type;
    static_assert(!sketch::wj::is_weighted_sketch<final_type>::value, "Can't have a weighted sketch be a final type");

    // Initialize vector of sketches
    std::vector<SketchType> sketches;
    sketches.reserve(inpaths.size());
    
    // Calculate sketch size and construct sketches
    const uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    while(sketches.size() < inpaths.size()) {
        sketches.emplace_back(construct<SketchType>(sketch_size));
        set_estim_and_jestim(sketches.back(), estim, jestim);
    }

    // Check if final sketch type is same as input sketch type
    static constexpr bool samesketch = std::is_same<SketchType, final_type>::value;
    final_type *final_sketches;
    std::unique_ptr<std::vector<final_type>> raii_final_sketches;

    std::atomic<uint32_t> ncomplete;
    ncomplete.store(0);
    const unsigned k = sp.k_;
    const unsigned wsz = sp.w_;
    RollingHasher<uint64_t> rolling_hasher(k, canon);
    if(inpaths.size() == 1 && presketched_only) {
        raii_final_sketches.reset(new std::vector<final_type>);
        gzFile ifp = gzopen(inpaths[0].data(), "rb");
        if(!ifp) UNRECOVERABLE_ERROR("Failed to open file.");
        for(;;) {
            raii_final_sketches->emplace_back(ifp);
            set_estim_and_jestim(raii_final_sketches->back(), estim, jestim);
        }
        final_sketches = raii_final_sketches->data();
        while(inpaths.size() < raii_final_sketches->size())
            inpaths.emplace_back(std::to_string(inpaths.size()));
        inpaths[0] = "0";
    } else {
        final_sketches =
            samesketch ? reinterpret_cast<final_type *>(sketches.data())
                       : static_cast<final_type *>(std::malloc(sizeof(*final_sketches) * inpaths.size()));
        OMP_PFOR_DYN
        for(size_t i = 0; i < sketches.size(); ++i) {
            const std::string &path(inpaths[i]);
            auto &sketch = sketches[i];
            if(presketched_only)  {
                CONST_IF(samesketch) {
                    sketch.read(path);
                    set_estim_and_jestim(sketch, estim, jestim); // HLL is the only type that needs this, and it's the same
                } else {
                    new(final_sketches + i) final_type(path.data()); // Read from path
                }
            } else {
                const std::string fpath(make_fname<SketchType>(path.data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix, enct));
                const bool isf = isfile(fpath);
                if(cache_sketch && isf) {
                    LOG_DEBUG("Sketch found at %s with size %zu, %u\n", fpath.data(), size_t(1ull << sketch_size), sketch_size);
                    CONST_IF(samesketch) {
                        sketch.read(fpath);
                        set_estim_and_jestim(sketch, estim, jestim);
                    } else {
                        new(final_sketches + i) final_type(fpath);
                    }
                } else {
                    const int tid = omp_get_thread_num();
                    Encoder<score::Lex> enc(nullptr, 0, sp, nullptr, canon);
                    auto ksp = &kseqs[tid];
                    if(cms.empty()) {
                        auto &h = sketch;
                        if(enct == BONSAI) for_each_substr([&](const char *s) {enc.for_each([&](u64 kmer){h.addh(kmer);}, s, ksp);}, inpaths[i], FNAME_SEP);
                        else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, ksp);}, inpaths[i], FNAME_SEP);
                        else for_each_substr([&](const char *s) {rolling_hasher.for_each_hash([&](u64 kmer){h.addh(kmer);}, s, ksp);}, inpaths[i], FNAME_SEP);
                    } else {
                        CountingSketch &cm = cms.at(tid);
                        const auto lfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) sketch.addh(kmer);};
                        if(enct == BONSAI)      for_each_substr([&](const char *s) {enc.for_each(lfunc, s, ksp);}, inpaths[i], FNAME_SEP);
                        else if(enct == NTHASH) for_each_substr([&](const char *s) {enc.for_each_hash(lfunc, s, ksp);}, inpaths[i], FNAME_SEP);
                        else                    for_each_substr([&](const char *s) {rolling_hasher.for_each_hash(lfunc, s, ksp);}, inpaths[i], FNAME_SEP);
                        cm.clear();
                    }
                    CONST_IF(!samesketch) new(final_sketches + i) final_type(std::move(sketch));
                    CONST_IF(samesketch) {
                        if(cache_sketch && !isf) sketch.write(fpath);
                    } else if(cache_sketch) final_sketches[i].write(fpath);
                }
            }
            ++ncomplete; // Atomic
        }
    }
    OMP_PFOR
    for(size_t i = 0; i < sketches.size(); ++i) {
        try {
            sketch_finalize(final_sketches[i]);
        } catch(const std::exception &ex) {
            std::cerr << "Failed to finalize sketch " << i << " for value " << inpaths[i] << ".\n"
            << "msg: " << ex.what() << '\n';
        }
    }
    kseqs.free();
    ks::string str("#Path\tSize (est.)\n");
    assert(str == "#Path\tSize (est.)\n");
    str.resize(BUFFER_FLUSH_SIZE);
    {
        const int fn(fileno(ofp));
        for(size_t i(0); i < sketches.size(); ++i) {
            double card;
            CONST_IF(samesketch) card = cardinality_estimate(sketches[i]);
            else                 card = cardinality_estimate(final_sketches[i]);
            str.sprintf("%s\t%zu\n", inpaths[i].data(), size_t(card));
            if(str.size() >= BUFFER_FLUSH_SIZE) str.flush(fn);
        }
        str.flush(fn);
    }
    if(ofp != stdout) std::fclose(ofp);
    str.clear();
    if(emit_fmt == UT_TSV && !nq) {
        str.sprintf("##Names\t");
        for(size_t i = 0; i < inpaths.size() - nq; ++i)
            str.sprintf("%s\t", inpaths[i].data());
        str.back() = '\n';
        str.write(fileno(pairofp)); str.free();
    } else if(emit_fmt == UPPER_TRIANGULAR) { // emit_fmt == UPPER_TRIANGULAR
        std::fprintf(pairofp, "%zu\n", inpaths.size());
        std::fflush(pairofp);
    }
    if(emit_fmt & NEAREST_NEIGHBOR_TABLE) {
        DBG_ONLY(std::fprintf(stderr, "[%s] About to make nn table with result type = %s. Number inpaths: %zu. \n", __PRETTY_FUNCTION__, emt2str(result_type), inpaths.size());)
        nndist_loop(pairofp, final_sketches, inpaths, k, result_type, emit_fmt, nq);
    } else {
        dist_loop<final_type>(pairofp, outpath, final_sketches, inpaths, use_scientific, k, result_type, emit_fmt, nthreads, BUFFER_FLUSH_SIZE, nq);
    }
    CONST_IF(!samesketch) {
        if(!raii_final_sketches) {
#if __cplusplus >= 201703L
            std::destroy_n(final_sketches, inpaths.size());
#else
            std::for_each(final_sketches, final_sketches + inpaths.size(), [](auto &sketch) {
                using T = typename std::decay<decltype(sketch)>::type;
                sketch.~T();
            });
#endif
            std::free(final_sketches);
        }
    }
} // dist_sketch_and_cmp


/**
 * @brief Macro to declare template specializations of dist_sketch_and_cmp for different sketch types
 * 
 * This macro declares three template specializations for the dist_sketch_and_cmp function:
 * 1. Base sketch type (DS)
 * 2. Weighted sketcher wrapping base type (WeightedSketcher<DS>)
 * 3. Weighted sketcher with exact counting adapter wrapping base type (WeightedSketcher<DS, ExactCountingAdapter>)
 *
 * Each specialization takes the same parameters:
 * @param inpaths Vector of input file paths
 * @param cms Vector of counting sketches
 * @param kseqs Buffer holder for k-mer sequences
 * @param ofp Output file pointer
 * @param pairofp Pair output file pointer reference
 * @param sp Spacer for k-mer sampling
 * @param ssarg Sketch size argument
 * @param mincount Minimum count threshold
 * @param estim Estimation method
 * @param jestim Joint estimation method
 * @param cache_sketch Whether to cache sketches
 * @param result_type Type of result to emit
 * @param emit_fmt Format to emit results in
 * @param presketched_only Whether to only use presketched data
 * @param nthreads Number of threads to use
 * @param use_scientific Whether to use scientific notation
 * @param suffix Suffix for output files
 * @param prefix Prefix for output files
 * @param canon Whether to canonicalize k-mers
 * @param entropy_minimization Whether to use entropy minimization
 * @param spacing Spacing pattern for k-mer sampling
 * @param nq Number of queries
 * @param enct Encoding type
 */
#define DECSKETCHCMP(DS) \
template void ::bns::dist_sketch_and_cmp<DS>(std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp,\
                   std::string,\
                   Spacer sp,\
                   unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                   bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                   size_t nq, EncodingType enct);\
template void ::bns::dist_sketch_and_cmp<sketch::wj::WeightedSketcher<DS>>(std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp,\
                   std::string,\
                   Spacer sp,\
                   unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                   bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                   size_t nq, EncodingType enct);\
template void ::bns::dist_sketch_and_cmp<sketch::wj::WeightedSketcher<DS, wj::ExactCountingAdapter>>(std::vector<std::string> &inpaths, std::vector<::bns::CountingSketch> &cms, KSeqBufferHolder &kseqs, std::FILE *ofp, std::FILE *&pairofp,\
                   std::string,\
                   Spacer sp,\
                   unsigned ssarg, unsigned mincount, EstimationMethod estim, JointEstimationMethod jestim, bool cache_sketch, EmissionType result_type, EmissionFormat emit_fmt,\
                   bool presketched_only, unsigned nthreads, bool use_scientific, std::string suffix, std::string prefix, bool canon, bool entropy_minimization, std::string spacing,\
                   size_t nq, EncodingType enct);

/**
 * @brief Flags controlling sketch behavior
 */
enum SketchFlags {
    SKIP_CACHED  = 1,  // Skip sketching if cached version exists
    CANONICALIZE = 2,  // Use canonical k-mers
    ENTROPY_MIN  = 4   // Use entropy minimization
};

/**
 * @brief Core sketching function template
 *
 * @tparam SketchType Type of sketch to use
 * @param ssarg Sketch size argument (log2 of bytes)
 * @param nthreads Number of threads to use
 * @param wsz Window size for k-mer sampling
 * @param k K-mer size
 * @param sp Spacer for k-mer sampling
 * @param inpaths Input sequence file paths
 * @param suffix Suffix for sketch filenames
 * @param prefix Prefix for sketch filenames
 * @param cms Vector of counting sketches
 * @param estim Method for cardinality estimation
 * @param jestim Method for joint cardinality estimation
 * @param kseqs Buffer for reading sequences
 * @param use_filter Vector indicating which sequences to filter
 * @param spacing K-mer sampling spacing pattern
 * @param sketchflags Flags controlling sketch behavior
 * @param mincount Minimum k-mer count threshold
 * @param enct Encoding type for k-mers
 * @param output_file Output file path for sketches
 */
template<typename SketchType>
INLINE void sketch_core(uint32_t ssarg, uint32_t nthreads, uint32_t wsz, uint32_t k, const Spacer &sp,
                        const std::vector<std::string> &inpaths, const std::string &suffix, const std::string &prefix,
                        std::vector<CountingSketch> &cms, EstimationMethod estim, JointEstimationMethod jestim,
                        KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,
                        int sketchflags, uint32_t mincount, EncodingType enct, std::string output_file)
{
    const auto canon = sketchflags & CANONICALIZE, skip_cached = sketchflags & SKIP_CACHED;
    std::vector<SketchType> sketches, sketchdest;
    uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    while(sketches.size() < (u32)nthreads) sketches.push_back(construct<SketchType>(sketch_size)), set_estim_and_jestim(sketches.back(), estim, jestim);
    if(output_file.size()) {
        sketchdest.reserve(inpaths.size());
        for(size_t i = 0, e = inpaths.size(); i < e; ++i) {
            sketchdest.push_back(construct<SketchType>(sketch_size));
            set_estim_and_jestim(sketchdest.back(), estim, jestim);
        }
    }
    std::vector<std::string> fnames(nthreads);
    RollingHasher<uint64_t> rolling_hasher(k, canon);

    gzFile outputfp = nullptr;
    if(output_file.size()) {
        if((outputfp = gzopen((output_file + ".labels.gz").data(), "w")) == nullptr) UNRECOVERABLE_ERROR("Failed to write sequence labels to file");
        for(const auto &path: inpaths) {
            gzwrite(outputfp, path.data(), path.size());
            gzputc(outputfp, '\n');
        }
        gzclose(outputfp);
        LOG_DEBUG("Wrote labels to file\n");
    }
    {
        if(outputfp) {
            if((outputfp = gzopen(output_file.data(), "w")) == nullptr)
                UNRECOVERABLE_ERROR("Failed to write sketches to file");
#ifndef NDEBUG
            else std::fprintf(stderr, "Opened file at %s\n", output_file.data());
#endif
        }
        OMP_PFOR_DYN
        for(size_t i = 0; i < inpaths.size(); ++i) {
            const int tid = omp_get_thread_num();
            std::string &fname = fnames[tid];
            fname = make_fname<SketchType>(inpaths[i].data(), sketch_size, wsz, k, sp.c_, spacing, suffix, prefix, enct);
            LOG_INFO("fname: %s from %s with tid = %d\n", fname.data(), inpaths[i].data(), tid);
            unsigned sketchind = outputfp ? unsigned(i): tid;
            auto &h = sketches[sketchind];
            if(skip_cached && isfile(fname)) {
                if(outputfp) h.read(fname);
                else continue;
            }
            Encoder<bns::score::Lex> enc(nullptr, 0, sp, nullptr, canon);
            const auto &path = inpaths[i];
            if(use_filter.size() && use_filter[i]) {
                auto &cm = cms[tid];
                auto lfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) h.add(kmer);};
                auto hlfunc = [&](u64 kmer){if(cm.addh(kmer) >= mincount) h.addh(kmer);};
// Macro to process each substring with a wrapper function and add kmers to sketch
#define FOR_EACH_FUNC(wrapperfunc) for_each_substr([&](const char *s) {wrapperfunc(lfunc, s, &kseqs[tid]);}, path, FNAME_SEP)

// Macro to process each substring with a wrapper function and add hashed kmers to sketch 
#define FOR_EACH_HASH_FUNC(wrapperfunc) for_each_substr([&](const char *s) {wrapperfunc(hlfunc, s, &kseqs[tid]);}, path, FNAME_SEP)

                // Process kmers based on encoding type
                switch(enct) {
                case NTHASH: FOR_EACH_FUNC(enc.for_each_hash); break;  // Use ntHash encoding
                case BONSAI: FOR_EACH_HASH_FUNC(enc.for_each); break;  // Use Bonsai encoding
                default: FOR_EACH_FUNC(rolling_hasher.for_each_hash); break; // Use default rolling hash
                }
                cm.clear(); // Clear counting sketch after use
            } else {
                // Define functions to add kmers directly to sketch without filtering
                auto lfunc = [&](u64 kmer){h.add(kmer);}; // Add kmer
                auto hlfunc = [&](u64 kmer){h.addh(kmer);}; // Add hashed kmer
                
                // Process kmers based on encoding type (same as above but without filtering)
                switch(enct) {
                case NTHASH: FOR_EACH_FUNC(enc.for_each_hash); break;
                case BONSAI: FOR_EACH_HASH_FUNC(enc.for_each); break;
                default: FOR_EACH_FUNC(rolling_hasher.for_each_hash); break;
                }
            }
#undef FOR_EACH_FUNC

            sketch_finalize(h); // Finalize the sketch after adding all kmers

            // Handle sketch output
            if(!outputfp) { 
                h.write(fname.data()); // Write sketch directly to disk
                h.clear(); // Clear sketch memory
            } else {        
                sketchdest[i] = h; // Copy to destination array
                h.reset(); // Reset sketch for reuse
            }
        }

        // Write all sketches to output file if specified
        if(outputfp) {
            DBG_ONLY(size_t j = 0; auto it = inpaths.cbegin();)
            for(const auto &sketch: sketchdest) {
                // Debug output showing progress
                DBG_ONLY(std::fprintf(stderr, "Writing sketch %zu/%zu from path %s to file\n", j++, inpaths.size(), (it++)->data());)
                sketch.write(outputfp); // Write sketch to output file
            }
            gzclose(outputfp); // Close output file
        }
    }
}

/**
 * @brief Core function for sketching sequences from a single input file
 *
 * @tparam SketchType Type of sketch to use (e.g. HLL, MinHash)
 * @param ssarg Sketch size argument (log2 of bytes)
 * @param nthreads Number of threads to use
 * @param sp Spacer for k-mer sampling
 * @param inpath Input sequence file path
 * @param outpath Output path for sketches
 * @param cm Counting sketch for filtering k-mers
 * @param estim Method for cardinality estimation
 * @param jestim Method for joint cardinality estimation
 * @param use_filter Whether to filter k-mers by count
 * @param sketchflags Flags controlling sketch behavior
 * @param mincount Minimum k-mer count threshold
 * @param enct Encoding type for k-mers
 */
template<typename SketchType>
INLINE void sketch_by_seq_core(uint32_t ssarg, uint32_t nthreads, const Spacer &sp,
                        const std::string &inpath, const std::string &outpath,
                        CountingSketch *cm, EstimationMethod estim, JointEstimationMethod jestim,
                        bool use_filter,
                        int sketchflags,
                        uint32_t mincount, EncodingType enct)
{
    // Extract canonicalization flag and calculate sketch size
    const auto canon = sketchflags & CANONICALIZE;
    const uint32_t sketch_size = bytesl2_to_arg(ssarg, SketchEnum<SketchType>::value);
    const auto k = sp.k_;

    // Initialize sketch and set estimation methods
    SketchType working_sketch(construct<SketchType>(sketch_size));
    set_estim_and_jestim(working_sketch, estim, jestim);

    // Initialize hashers and encoder
    RollingHasher<uint64_t> rolling_hasher(k, canon);
    Encoder<bns::score::Lex> enc(nullptr, 0, sp, nullptr, canon);

    // Open input file
    gzFile fp = gzopen(inpath.data(), "rb");
    if(!fp) throw ZlibError(std::string("Failed to open file for reading at ") + inpath);

    // Open output files for sketches and sequence names
    gzFile ofp = gzopen(outpath.data(), "wb");
    if(!ofp) throw ZlibError(std::string("Failed to open file for writing at ") + outpath);

    std::string namepath = outpath == "/dev/stdout"
        ? std::string("stdout.names")
        : outpath + ".names";
    std::FILE *nameofp = fopen(namepath.data(), "w");
    if(!nameofp) UNRECOVERABLE_ERROR(std::string("Failed to open file for writing at ") + namepath);
    fprintf(nameofp, "#k=%d:Names for sequences sketched\n", k);

    // Initialize sequence reader
    kseq_t *ks = kseq_init(fp);
    auto &h = working_sketch;

    // Define k-mer processing functions
    auto add = [&](u64 kmer) {h.addh(kmer);};
    auto cadd = [&](u64 kmer){if(cm->addh(kmer) >= mincount) h.addh(kmer);};

    // Process each sequence
    std::vector<std::string> seqnames;
    while(kseq_read(ks) >= 0) {
        // Process k-mers with or without filtering
        if(use_filter) {
            if(enct == NTHASH)
                enc.for_each_hash(cadd, ks->seq.s, ks->seq.l);
            else if(enct == BONSAI)
                enc.for_each(cadd, ks->seq.s, ks->seq.l);
            else
                rolling_hasher.for_each_hash(cadd, ks->seq.s, ks->seq.l);
            cm->clear();
        } else {
            if(enct == NTHASH)
                enc.for_each_hash(add, ks->seq.s, ks->seq.l);
            else if(enct == BONSAI)
                enc.for_each(add, ks->seq.s, ks->seq.l);
            else
                rolling_hasher.for_each_hash(add, ks->seq.s, ks->seq.l);
        }

        // Finalize sketch and handle errors
        try {
            sketch_finalize(h);
        } catch(const std::exception &ex) {
            std::cerr << "Failed to finalize sketch  for sequence " << ks->seq.s << ".\n"
            << "msg: " << ex.what() << '\n';
        }

        // Write sequence name and sketch
        if(std::fwrite(ks->name.s, 1, ks->name.l, nameofp) != ks->name.l) std::fprintf(stderr, "Warning: error in writing sequence name\n");
        std::fputc('\n', nameofp);
        h.write(ofp);
        h.clear();
    }

    // Clean up
    kseq_destroy(ks);
    gzclose(fp);
    gzclose(ofp);
    std::fclose(nameofp);
}

// Type alias for storing a distance/similarity value and corresponding index
using validx_t = std::pair<float, uint32_t>;

/**
 * @brief Updates nearest neighbor data with thread synchronization
 * 
 * @tparam Cmp Type of comparison function
 * @param val Distance/similarity value to potentially insert
 * @param ptr Array of current nearest neighbors
 * @param nneighbors Number of neighbors to track
 * @param j Index of current comparison
 * @param mut Mutex for thread synchronization
 * @param cmp Comparison function for ordering neighbors
 */
template<typename Cmp>
INLINE void lock_update(float val, validx_t *ptr,
                        unsigned nneighbors,
                        size_t j, std::mutex &mut, const Cmp &cmp)
{
    // Check if new value should be inserted based on comparison
    if(cmp(val, ptr->first)) {
        // Lock to prevent concurrent updates
        std::lock_guard<std::mutex> lg(mut);
        // Re-check after acquiring lock in case another thread updated
        if(cmp(val, ptr->first)) {
            // Remove worst neighbor
            std::pop_heap(ptr, ptr + nneighbors, cmp);
            // Insert new neighbor
            ptr[nneighbors - 1] = {val, j};
            // Restore heap property
            std::push_heap(ptr, ptr + nneighbors, cmp);
        }
    }
}

/**
 * @brief Updates nearest neighbor data without thread synchronization
 *
 * @param val Distance/similarity value to potentially insert
 * @param ptr Array of current nearest neighbors
 * @param nneighbors Number of neighbors to track
 * @param j Index of current comparison
 * @param is_greater True if using similarity measure (greater is better), false for distance
 */
INLINE void lockfree_update(float val, validx_t *ptr,
                            unsigned nneighbors,
                            size_t j, bool is_greater)
{
    if(is_greater) {
        // For similarity measures, keep larger values
        if(val > ptr->first) {
            std::pop_heap(ptr, ptr + nneighbors, std::greater<>());
            ptr[nneighbors - 1] = validx_t(val, j);
            std::push_heap(ptr, ptr + nneighbors, std::greater<>());
        }
    } else {
        // For distance measures, keep smaller values
        if(val < ptr->first) {
            std::pop_heap(ptr, ptr + nneighbors, std::less<>());
            ptr[nneighbors - 1] = validx_t(val, j);
            std::push_heap(ptr, ptr + nneighbors, std::less<>());
        }
    }
}

/**
 * @brief Updates nearest neighbor data with thread synchronization
 *
 * @param val Distance/similarity value to potentially insert
 * @param ptr Array of current nearest neighbors
 * @param nneighbors Number of neighbors to track
 * @param j Index of current comparison
 * @param mut Mutex for thread synchronization
 * @param cmp Comparison function for ordering neighbors
 */
template<typename Cmp>
INLINE void lock_update(float val, validx_t *ptr,
                        unsigned nneighbors,
                        size_t j, std::mutex &mut, const Cmp &cmp)
{
    if(cmp(val, ptr->first)) {
        std::lock_guard<std::mutex> lg(mut);
        if(cmp(val, ptr->first)) { // after getting the lock, check again
            std::pop_heap (ptr, ptr + nneighbors, cmp);
            ptr[nneighbors - 1] = {val, j};
            std::push_heap(ptr, ptr + nneighbors, cmp);
        }
    }
}

/**
 * @brief Updates nearest neighbor data without thread synchronization
 *
 * @param val Distance/similarity value to potentially insert
 * @param ptr Array of current nearest neighbors
 * @param nneighbors Number of neighbors to track
 * @param j Index of current comparison
 * @param is_greater True if using similarity measure (greater is better), false for distance
 */
INLINE void lockfree_update(float val, validx_t *ptr,
                            unsigned nneighbors,
                            size_t j, bool is_greater)
{
    if(is_greater) {
        if(val > ptr->first) {
            std::pop_heap (ptr, ptr + nneighbors, std::greater<>());
            ptr[nneighbors - 1] = validx_t(val, j);
            std::push_heap(ptr, ptr + nneighbors, std::greater<>());
        }
    } else {
        if(val < ptr->first) {
            std::pop_heap (ptr, ptr + nneighbors, std::less<>());
            ptr[nneighbors - 1] = validx_t(val, j);
            std::push_heap(ptr, ptr + nneighbors, std::less<>());
        }
    }
}

/**
 * @brief Performs nearest neighbor search between sketches
 *
 * @param neighbors Output array to store nearest neighbor results
 * @param sketches Array of sketch data structures
 * @param inpaths Vector of input file paths
 * @param k K-mer size used for sketching
 * @param result_type Type of distance/similarity measure
 * @param nq Number of query sequences (0 for all-vs-all comparison)
 * @param nneighbors Number of nearest neighbors to find per sequence
 * @param is_similarity True if using similarity measure, false for distance
 * @param func Function to compute distance/similarity between sketches
 */
template<typename SketchType, typename Func>
void perform_nns(validx_t *neighbors,
                 SketchType *sketches, const std::vector<std::string> &inpaths,
                 const unsigned k, const EmissionType result_type,
                 size_t nq,
                 unsigned nneighbors,
                 const bool is_similarity, const Func &func)
{
    // Calculate number of sequences to process
    const size_t n = nq ? nq: inpaths.size();

    // Initialize neighbor arrays with extreme values
    {
        const validx_t empty_val(is_similarity ? -std::numeric_limits<float>::max(): std::numeric_limits<float>::max(),
                                 uint32_t(-1));
        OMP_PFOR
        for(size_t i = 0; i < n; ++i) {
            std::fill_n(&neighbors[i * nneighbors], nneighbors,
                        empty_val);
        }
    }

    if(nq == 0) {
        // All-vs-all comparison case
        auto mutexes = std::make_unique<std::mutex[]>(n);
        OMP_PFOR_DYN
        for(size_t i = 0; i < n; ++i) {
            auto lhptr = &neighbors[i * nneighbors];
            const auto &h1 = sketches[i];
            for(size_t j = i + 1; j < n; ++j) {
                auto rhptr = &neighbors[j * nneighbors];
                const float val = func(sketches[j], h1);
                if(is_similarity) {
                    lock_update(val, lhptr, nneighbors, j, mutexes[i], std::greater<>());
                    lock_update(val, rhptr, nneighbors, i, mutexes[j], std::greater<>());
                } else {
                    lock_update(val, lhptr, nneighbors, j, mutexes[i], std::less<>());
                    lock_update(val, rhptr, nneighbors, i, mutexes[j], std::less<>());
                }
            }
        }
    } else {
        // Query-vs-reference comparison case
        const size_t npaths = inpaths.size(), nr = npaths - nq;
        OMP_PFOR_DYN
        for(size_t qi = nr; qi < npaths; ++qi) {
            const size_t qind = qi - nr;
            const auto srcptr = &neighbors[qind * nneighbors];
            const auto &h1 = sketches[qi];
            for(size_t j = 0; j < nr; ++j) {
                lockfree_update(func(sketches[j], h1), srcptr, nneighbors, j, is_similarity);
            }
        }
    }

    LOG_DEBUG("Finished loop, now sorting\n");
    
    // Sort final neighbor lists
    OMP_PFOR
    for(size_t i = 0; i < n; ++i) {
        auto start = neighbors + (i * nneighbors), end = start + nneighbors;
        is_similarity ? std::sort(start, end, std::greater<>())
                      : std::sort(start, end, std::less<>());
    }
}

/**
 * @brief Performs core distance/similarity computation between sketches
 *
 * @tparam MT Whether to use multithreading (default: true)
 * @tparam SketchType Type of sketch being compared
 * @tparam T Type of distance/similarity values
 * @tparam Func Type of comparison function
 * @param dists Array to store computed distances
 * @param nsketches Number of sketches to compare
 * @param sketches Array of sketches
 * @param func Comparison function to use
 * @param i Index of current sketch being compared
 */
template<bool MT=true, typename SketchType, typename T, typename Func>
inline void perform_core_op(T *dists, size_t nsketches, SketchType *sketches, const Func &func, size_t i) {
    auto &h1 = sketches[i];
#define compute_j(j) do {dists[j - i - 1] = func(sketches[j], h1);} while(0)
    if(MT) {
        OMP_PFOR_DYN
        for(size_t j = i + 1; j < nsketches; ++j) compute_j(j);
    } else {
        for(size_t j = i + 1; j < nsketches; ++j) compute_j(j);
    }
#undef compute_j
}

/**
 * @brief Computes and outputs nearest neighbor distances between sketches
 *
 * @tparam SketchType Type of sketch being compared
 * @param ofp Output file pointer
 * @param sketches Array of sketches
 * @param inpaths Vector of input paths
 * @param k K-mer size used
 * @param result_type Type of distance/similarity measure
 * @param emit_fmt Output format
 * @param nq Number of query sketches (0 for all-vs-all)
 */
template<typename SketchType>
void nndist_loop(std::FILE *ofp, SketchType *sketches,
               const std::vector<std::string> &inpaths,
               const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt,
               size_t nq) {
    // Display global arguments
    gargs.show();

    // Get number of neighbors to find
    unsigned nneighbors = gargs.number_neighbors;
    size_t npairs = (nq ? nq: inpaths.size());
    size_t possible_num_neighbors = nq ? inpaths.size() - nq: inpaths.size();

    // Adjust number of neighbors if more requested than possible
    if(nneighbors > possible_num_neighbors) {
        std::fprintf(stderr, "Only reporting %zu rather than %u neighbors due to their being only that many sets.\n", possible_num_neighbors, nneighbors);
        nneighbors = possible_num_neighbors;
    }

    // Calculate offsets and allocate neighbor array
    const size_t qoffset = nq ? inpaths.size() - nq: size_t(0);
    const size_t ntups = npairs * nneighbors;
    auto neighbors = std::make_unique<validx_t[]>(ntups);
    std::fprintf(stderr, "made %zu pairs with %zu nq, %u neighbors and %zu tups\n", npairs, nq, nneighbors, ntups);

    // Set up comparison function
    const double ksinv = 1./ k;
    auto call_cmp = [result_type, ksinv](const auto &x, const auto &y) {return result_cmp(x, y, result_type, ksinv);};

    // Perform nearest neighbor search
    perform_nns(neighbors.get(), sketches, inpaths, k, result_type, nq, nneighbors, emt2nntype(result_type) == SIMILARITY_MEASURE, call_cmp);

    // Output results in binary format
    if(emit_fmt & BINARY) {
        uint32_t n = inpaths.size();
        std::fwrite(&n, sizeof(n), 1, ofp);
        n = nneighbors;
        std::fwrite(&n, sizeof(n), 1, ofp);
        size_t nb = nneighbors * inpaths.size();
        if(unlikely(std::fwrite(neighbors.get(), sizeof(validx_t), nb, ofp) != nb))
            UNRECOVERABLE_ERROR("Failed to write neighbors to disk (binary)\n");
    }
    // Output results in text format 
    else {
        std::fprintf(ofp, "#File\tNeighbor ID:distance\t...\n");
        
        // Get number of threads
        int nt = 1;
#ifdef _OPENMP
        _Pragma("omp parallel")
        {
            _Pragma("omp single")
            nt = omp_get_num_threads();
        }
#endif
        // Create string buffers for each thread
        std::vector<ks::string> kstrs;
        while(kstrs.size() < unsigned(nt))
            kstrs.emplace_back(1ull<<10);
        std::fprintf(stderr, "Made buffer, ones per thread\n");

        const size_t npaths = inpaths.size();
        
        // Write neighbor data in parallel
        OMP_PFOR
        for(size_t i = 0; i < npaths; ++i) {
            auto tid = OMP_ELSE(omp_get_thread_num(), 0);
            auto &buf(kstrs[tid]);
            const validx_t *nptr = &neighbors[i * nneighbors];
            assert(i * nneighbors + nneighbors < npairs);
            size_t nameind = i + qoffset;
            assert(nameind < inpaths.size());
            buf += inpaths[nameind];
            
            // Write each neighbor's ID and distance
            for(unsigned j = 0; j < nneighbors; ++j) {
                assert(nptr + j < neighbors.get() + ntups || std::fprintf(stderr, "i = %zu, j = %zu, offset = %zu\n", i, size_t(j), i * nneighbors + j));
                const auto ndat = nptr[j];
                buf.putc_('\t');
                buf.putw_(ndat.second);
                buf.putc_(':');
                buf.sprintf("%g", ndat.first);
            }
            buf.putc_('\n');
            
            // Flush buffer if it gets too large
            if(buf.size() > (1u << 14)) {
                OMP_CRITICAL
                {
                    buf.flush(ofp);
                }
            }
        }
        
        // Final flush of all buffers
        for(auto buf: kstrs) buf.flush(ofp);
    }
}

/**
 * @brief Main distance calculation loop between sketches
 *
 * @param ofp Output file pointer for writing distances
 * @param ofp_name Name of output file 
 * @param sketches Array of sketch data structures
 * @param inpaths Vector of input file paths
 * @param use_scientific Whether to use scientific notation in output
 * @param k K-mer size used for sketching
 * @param result_type Type of distance/similarity metric to compute
 * @param emit_fmt Format for emitting results (e.g. TSV, binary)
 * @param unused Unused parameter
 * @param buffer_flush_size Size at which to flush output buffer
 * @param nq Number of query sequences (0 for all-vs-all comparison)
 */
template<typename SketchType>
void dist_loop(std::FILE *&ofp, std::string ofp_name, SketchType *sketches, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, int, const size_t buffer_flush_size, size_t nq)
{
    // Handle query vs reference case
    if(nq) {
        partdist_loop<SketchType>(ofp, sketches, inpaths, use_scientific, k, result_type, emit_fmt, buffer_flush_size, nq);
        return;
    }

    // Validate symmetric distance calculation
    if(!is_symmetric(result_type)) {
        char buf[1024];
        std::sprintf(buf, "Can't perform symmetric distance comparisons with a symmetric method (%s/%d). To perform an asymmetric distance comparison between a given set and itself, provide the same list of filenames to both -Q and -F.\n", emt2str(result_type), int(result_type));
        UNRECOVERABLE_ERROR(buf);
    }

    // Initialize parameters
    const float ksinv = 1./ k;
    const int pairfi = fileno(ofp);
    const size_t nsketches = inpaths.size();
    std::future<size_t> submitter;
    
    // Define comparison function
    auto cmp = [ksinv,result_type](const auto &x, const auto &y) {return result_cmp(x, y, result_type, ksinv);};

    // Handle text output formats
    if((emit_fmt & BINARY) == 0) {
        // Allocate buffers for distances
        std::array<std::vector<float>, 2> dps;
        dps[0].resize(nsketches - 1);
        dps[1].resize(std::max(ssize_t(nsketches) - 2, ssize_t(1)));
        ks::string str;

        // Calculate and emit distances for each sketch
        for(size_t i = 0; i < nsketches; ++i) {
            auto distp = dps[i & 1].data();
            perform_core_op<true>(distp, nsketches, sketches, cmp, i);
            if(i) submitter.get();
            submitter = std::async(std::launch::async, submit_emit_dists<>,
                                   pairfi, distp, nsketches, i,
                                   std::ref(str), std::ref(inpaths), emit_fmt, use_scientific, buffer_flush_size);
        }
        submitter.get();
    }
    // Handle binary output formats 
    else {
        const float defv = static_cast<float>(emt2nntype(result_type) == SIMILARITY_MEASURE);
        
        // Binary matrix format
        if(emit_fmt == BINARY) {
            // Handle terminal output
            if(::isatty(::fileno(ofp))) {
                std::thread sub;
                uint64_t ns = nsketches;
                std::fputc('\0', ofp), std::fwrite(&ns, sizeof(ns), 1, ofp);
                
                // Calculate distances row by row
                for(size_t i = 0; i < nsketches - 1; ++i) {
                    std::unique_ptr<float[]> dat(new float[nsketches - i - 1]);
                    OMP_PFOR
                    for(size_t j = i + 1; j < nsketches; ++j)
                        dat[j - i - 1] = cmp(sketches[i], sketches[j]);
                    if(sub.joinable()) sub.join();
                    sub = std::thread([dat=std::move(dat),n=nsketches - i - 1,ofp,i]() {
                        if(std::fwrite(dat.get(), sizeof(float), n, ofp) != n)
                            throw std::runtime_error(std::string("Failed to write row ") + std::to_string(i) + " to disk");
                    });
                }
                sub.join();
            }
            // Handle file output
            else {
                // Initialize file header
                std::rewind(ofp);
                std::fputc('\0', ofp);
                uint64_t nelem = nsketches;
                if(std::fwrite(&nelem, sizeof(nelem), 1, ofp) != 1) throw std::runtime_error("Failure");
                
                // Set file size
                const size_t fsz = 1 + sizeof(uint64_t) + ((nsketches * (nsketches - 1)) >> 1) * sizeof(float);
                ::ftruncate(::fileno(ofp), fsz);
                std::fclose(ofp);
                ofp = nullptr;
                
                // Fill distance matrix in parallel
                dm::DistanceMatrix<float> dm(ofp_name.data(), nsketches, defv);
                dm::parallel_fill(dm, nsketches, [&cmp,sketches](size_t i, size_t j) {return cmp(sketches[i], sketches[j]);}, gargs.nperbatch);
            }
        }
        // Full TSV format
        else {
            if(emit_fmt != FULL_TSV) throw std::runtime_error("Invalid emit_fmt");
            
            // Write header row
            std::fputs("#Names", ofp);
            for(size_t i = 0; i < inpaths.size(); ++i) {
                std::fputs(inpaths[i].data(), ofp);
                std::fputc(i == inpaths.size() - 1 ? '\n': '\t', ofp);
            }

            // Calculate and write distances row by row
            std::thread sub;
            for(size_t i = 0; i < nsketches; ++i) {
                std::unique_ptr<float[]> dat(new float[nsketches]);
                OMP_PFOR_DYN
                for(size_t j = 0; j < nsketches; ++j) {
                    if(j == i) dat[j] = 0.;
                    else       dat[j] = cmp(sketches[i], sketches[j]);
                }
                if(sub.joinable()) sub.join();
                sub = std::thread([dat=std::move(dat),nsketches,ofp,i,name=inpaths[i].data()]() {
                    std::fprintf(ofp, "%s\t", name);
                    size_t j;
                    for(j = 0; j < nsketches - 1; ++j) {
                        std::fprintf(ofp, "%0.6g\t", dat[j]);
                    }
                    std::fprintf(ofp, "%0.6g\n", dat[j]);
                });
            }
            if(sub.joinable()) sub.join();
        }
    }
} // dist_loop
/**
 * @brief Macro to declare template specializations of sketch_core for different sketch types
 * 
 * @param DS Base sketch type to generate specializations for
 *
 * Declares three template specializations:
 * 1. Base sketch type
 * 2. Weighted sketcher wrapping base type 
 * 3. Weighted sketcher with exact counting adapter wrapping base type
 *
 * Each specialization takes the following parameters:
 * @param ssarg Sketch size argument (log2 of bytes)
 * @param nthreads Number of threads to use
 * @param wsz Window size for k-mer sampling
 * @param k K-mer size
 * @param sp Spacer for k-mer sampling
 * @param inpaths Input sequence file paths
 * @param suffix Suffix for sketch filenames
 * @param prefix Prefix for sketch filenames
 * @param counting_sketches Vector of counting sketches
 * @param estim Method for cardinality estimation
 * @param jestim Method for joint cardinality estimation
 * @param kseqs Buffer for reading sequences
 * @param use_filter Vector indicating which sequences to filter
 * @param spacing K-mer sampling spacing pattern
 * @param sketchflags Flags controlling sketch behavior
 * @param mincount Minimum k-mer count threshold
 * @param enct Encoding type for k-mers
 * @param s Additional string parameter
 */
#define DECSKETCHCORE(DS) \
  template void sketch_core<DS>(uint32_t ssarg, uint32_t nthreads,\
                                uint32_t wsz, uint32_t k, const Spacer &sp,\
                                const std::vector<std::string> &inpaths,\
                                const std::string &suffix,\
                                const std::string &prefix, std::vector<CountingSketch> &counting_sketches,\
                                EstimationMethod estim, JointEstimationMethod jestim,\
                                KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,\
                                int sketchflags, uint32_t mincount, EncodingType enct, std::string s);\
  template void sketch_core<wj::WeightedSketcher<DS>>(uint32_t ssarg, uint32_t nthreads,\
                                uint32_t wsz, uint32_t k, const Spacer &sp,\
                                const std::vector<std::string> &inpaths,\
                                const std::string &suffix,\
                                const std::string &prefix, std::vector<CountingSketch> &counting_sketches,\
                                EstimationMethod estim, JointEstimationMethod jestim,\
                                KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,\
                                int sketchflags, uint32_t mincount, EncodingType enct, std::string s);\
  template void sketch_core<wj::WeightedSketcher<DS, wj::ExactCountingAdapter>>(uint32_t ssarg, uint32_t nthreads,\
                                uint32_t wsz, uint32_t k, const Spacer &sp,\
                                const std::vector<std::string> &inpaths,\
                                const std::string &suffix,\
                                const std::string &prefix, std::vector<CountingSketch> &counting_sketches,\
                                EstimationMethod estim, JointEstimationMethod jestim,\
                                KSeqBufferHolder &kseqs, const std::vector<bool> &use_filter, const std::string &spacing,\
                                int sketchflags, uint32_t mincount, EncodingType enct, std::string s);\


} // namespace bns
