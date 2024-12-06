#ifndef DASHING_H__
#define DASHING_H__
#include <omp.h>
#include "sketch/bbmh.h"
#include "sketch/mh.h"
#include "sketch/hmh.h"
#include "sketch/mult.h"
#include "sketch/hk.h"
#include "sketch/bf.h"
#include "bonsai/encoder.h"
#include "khset/khset.h"
#include "distmat/distmat.h"
#include <sstream>
#include "getopt.h"
#include <sys/stat.h>
#include "substrs.h"
#include "khset64.h"
#include "enums.h"

#ifndef _OPENMP
#error("Need OpenMP")
#endif

#if __cplusplus >= 201703L && __cpp_lib_execution
#include <execution>
#endif
#ifndef BUFFER_FLUSH_SIZE
#define BUFFER_FLUSH_SIZE (1u << 18)
#endif

#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},

#define SHARED_OPTS \
    LO_FLAG("wj-exact", 145, gargs.exact_weighted, true)\
    LO_FLAG("use-wide-hll", 144, sketch_type, WIDE_HLL) \
    LO_FLAG("defer-hll", 146, gargs.defer_hll_creation, true)\
    /*LO_FLAG("use-hyperminhash", 147, sketch_type, HYPERMINHASH)*/\


using BKHash64 = sketch::minhash::BottomKHasher<sketch::WangHash, uint64_t>;

#define DIST_LONG_OPTS \
static option_struct dist_long_options[] = {\
    LO_FLAG("avoid-sorting", 'n', avoid_fsorting, true)\
    LO_FLAG("by-entropy", 'g', entropy_minimization, true) \
    LO_FLAG("cache-sketches", 'W', cache_sketch, true)\
    LO_FLAG("countmin", 'y', sm, CBF)\
    LO_FLAG("emit-binary", 'b', emit_fmt, BINARY)\
    LO_FLAG("full-mash-dist", 'l', result_type, FULL_MASH_DIST)\
    LO_FLAG("full-tsv", 'T', emit_fmt, FULL_TSV)\
    LO_FLAG("no-canon", 'C', canon, false)\
    LO_FLAG("phylip", 'U', emit_fmt, UPPER_TRIANGULAR)\
    LO_FLAG("presketched", 'H', presketched_only, true)\
    LO_FLAG("sizes", 'Z', result_type, SIZES)\
    LO_FLAG("sketch-by-fname", 'N', sm, BY_FNAME)\
    LO_FLAG("use-bb-minhash", '8', sketch_type, BB_MINHASH)\
    LO_FLAG("use-scientific", 'e', use_scientific, true)\
    LO_ARG("bbits", 'B')\
    LO_ARG("cm-sketch-size", 't')\
    LO_ARG("ertl-joint-mle", 'J')\
    LO_ARG("ertl-mle", 'm')\
    LO_ARG("improved", 'I')\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("min-count", 'c')\
    LO_ARG("nhashes", 'q')\
    LO_ARG("nthreads", 'p')\
    LO_ARG("original", 'E')\
    LO_ARG("out-dists", 'O') \
    LO_ARG("out-sizes", 'o') \
    LO_ARG("paths", 'F')\
    LO_ARG("prefix", 'P')\
    LO_ARG("query-paths", 'Q') \
    LO_ARG("seed", 'R')\
    LO_ARG("sketch-size", 'S')\
    LO_ARG("spacing", 's')\
    LO_ARG("suffix", 'x')\
    LO_ARG("window-size", 'w')\
    LO_ARG("help", 'h')\
    /*LO_ARG("mkdist", 1337)*/\
    LO_FLAG("use-range-minhash", 128, sketch_type, RANGE_MINHASH)\
    LO_FLAG("use-full-khash-sets", 130, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("use-full-hash-sets", 1000, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("use-hash-sets", 1000, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("hash-sets", 10001, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("use-full-sets", 1000, sketch_type, FULL_KHASH_SET)\
    LO_FLAG("full-containment-dist", 133, result_type, FULL_CONTAINMENT_DIST) \
    LO_FLAG("use-bloom-filter", 134, sketch_type, BLOOM_FILTER)\
    LO_FLAG("use-nthash", 136, enct, NTHASH)\
    LO_FLAG("containment-index", 131, result_type, CONTAINMENT_INDEX) \
    LO_FLAG("containment-dist", 132, result_type, CONTAINMENT_DIST) \
    LO_FLAG("mash-dist", 'M', result_type, MASH_DIST)\
    LO_FLAG("symmetric-containment-index", 137, result_type, SYMMETRIC_CONTAINMENT_INDEX) \
    LO_FLAG("symmetric-containment-dist", 138, result_type, SYMMETRIC_CONTAINMENT_DIST) \
    LO_FLAG("use-cyclic-hash", 139, enct, CYCLIC)\
    LO_ARG("wj-cm-sketch-size", 140)\
    LO_ARG("wj-cm-nhashes", 141)\
    LO_FLAG("wj", 142, weighted_jaccard, true)\
    LO_ARG("nearest-neighbors", 143)\
    SHARED_OPTS \
    LO_ARG("nperbatch", 148)\
    {0,0,0,0}\
};


namespace bns {
// Forward declaration for khash set type
struct khset64_t;

// Import HyperLogLog and HyperMinHash types from sketch namespace
using sketch::mh::HyperLogLogHasher;
using sketch::HyperMinHash;

// Function to flatten multiple files into one output file
int flatten_all(const std::vector<std::string> &fpaths, const std::string outpath, std::vector<unsigned> &k_values);

namespace detail {
    // Helper function to sort paths by file size
    void sort_paths_by_fsize(std::vector<std::string> &paths);
}

// Get total size of files in a path
size_t posix_fsizes(const std::string &path, const char sep=FNAME_SEP);

// Import commonly used namespaces
using namespace sketch;
using namespace hll;

// Import specific types from sketch namespace
using sketch::BBitMinHasher;
using option_struct = struct option;
using sketch::WangHash;

// Type aliases for specialized MinHash variants
using CRMFinal = mh::FinalCRMinHash<uint64_t, uint32_t>;
using RMFinal = mh::FinalRMinHash<uint64_t, sketch::common::Allocator<uint64_t>>;

// Hash functor that applies a seed to the input before hashing
template<typename BaseHash>
struct SeededHash {
    BaseHash wh_;
    const uint64_t seed_;
    SeededHash(uint64_t seed): seed_(seed) {}
    uint64_t operator()(uint64_t x) const {return wh_(x ^ seed_);}
};


#if DASHING_USE_HK
#define DASHING_COUNTING_SKETCH ::sketch::hk::HeavyKeeper<6, 10, SeededHash<sketch::common::WangHash>>
#else
#define DASHING_COUNTING_SKETCH ::sketch::ccm_t
#endif
using CountingSketch = DASHING_COUNTING_SKETCH;
// Generic similarity function that returns Jaccard index between two sketches
template<typename T> inline double similarity(const T &a, const T &b) {
    return a.jaccard_index(b);
}

// Specialized similarity function for CRMFinal type that uses histogram intersection
template<> inline double similarity<CRMFinal>(const CRMFinal &a, const CRMFinal &b) {
    return a.histogram_intersection(b);
}

// Generic sketch finalization function that does nothing by default
template<typename T>
inline void sketch_finalize(T &x) {}

// Specialized sketch finalization for khset64_t type that converts to SHS format
template<> inline void sketch_finalize<khset64_t>(khset64_t &x) {x.cvt2shs();}

// Calculate distance index between two floating point values
// ji: Jaccard index
// ksinv: Inverse of k-mer size
// Returns: Distance metric based on Jaccard index
template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
// Calculate distance metric based on Jaccard index using log approximation
// ji: Jaccard index between two sketches
// ksinv: Inverse of k-mer size (1/k)
// Returns: Distance metric in range [0,1], where 0 means identical and 1 means completely different
typename std::common_type<FType1, FType2>::type dist_index(FType1 ji, FType2 ksinv) {
    return ji ? -std::log(2. * ji / (1. + ji)) * ksinv: 1.;
}

// Calculate distance metric based on containment using log approximation
// containment: Containment index between two sketches (|A ∩ B|/|A|)
// ksinv: Inverse of k-mer size (1/k)
// Returns: Distance metric in range [0,1], where 0 means complete containment and 1 means no containment
template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
typename std::common_type<FType1, FType2>::type containment_dist(FType1 containment, FType2 ksinv) {
    return containment ? -std::log(containment) * ksinv: 1.;
}

// Calculate distance metric based on Jaccard index without log approximation
// ji: Jaccard index between two sketches
// ksinv: Inverse of k-mer size (1/k)
// Returns: Distance metric in range [0,1] using direct power calculation
template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
// Calculate exact distance metric based on Jaccard index without log approximation
// ji: Jaccard index between two sketches
// ksinv: Inverse of k-mer size (1/k)
// Returns: Distance metric in range [0,1] using direct power calculation
typename std::common_type<FType1, FType2>::type full_dist_index(FType1 ji, FType2 ksinv) {
    return 1. - std::pow(2.*ji/(1. + ji), ksinv);
}

// Calculate exact distance metric based on containment without log approximation
// containment: Containment index between two sketches (|A ∩ B|/|A|)
// ksinv: Inverse of k-mer size (1/k)
// Returns: Distance metric in range [0,1] using direct power calculation
template<typename FType1, typename FType2,
         typename=typename std::enable_if<
            std::is_floating_point<FType1>::value && std::is_floating_point<FType2>::value
          >::type
         >
typename std::common_type<FType1, FType2>::type full_containment_dist(FType1 containment, FType2 ksinv) {
    return 1. - std::pow(containment, ksinv);
}

// Calculate containment index between two sketches
// Returns: |A ∩ B|/|A| where A and B are the sets represented by the sketches
template<typename T>
inline double containment_index(const T &a, const T &b) {
    return a.containment_index(b);
}

// Calculate full set comparison metrics between two sketches
// Returns: Array containing intersection size, size of first set, size of second set
template<typename T>
inline std::array<double, 3> set_triple(const T &a, const T &b) {
    return a.full_set_comparison(b);
}

// Prevent calling set_triple on weighted sketches before finalization
template<typename T>
inline std::array<double, 3> set_triple(const wj::WeightedSketcher<T> &a, const wj::WeightedSketcher<T> &b) {
    UNRECOVERABLE_ERROR("This should only be called on the finalized sketches");
}

// Macro to define error cases for sketch types that don't support containment operations
#define CONTAIN_OVERLOAD_FAIL(x)\
template<>\
inline double containment_index<x>(const x &b, const x &a) {\
    UNRECOVERABLE_ERROR(std::string("Containment index not implemented for ") + __PRETTY_FUNCTION__);\
}\
template<>\
inline std::array<double, 3> set_triple<x>(const x &b, const x &a) {\
    UNRECOVERABLE_ERROR(std::string("set_triple not implemented for ") + __PRETTY_FUNCTION__);\
}

// Define error cases for specific sketch types
CONTAIN_OVERLOAD_FAIL(RMFinal)
CONTAIN_OVERLOAD_FAIL(bf::bf_t)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<RMFinal>)
CONTAIN_OVERLOAD_FAIL(wj::WeightedSketcher<bf::bf_t>)
CONTAIN_OVERLOAD_FAIL(CRMFinal)
// Define weighted sketch types using exact counting adapter
using wjRMFinal = wj::WeightedSketcher<RMFinal, wj::ExactCountingAdapter>;
using wjBFFinal = wj::WeightedSketcher<bf::bf_t, wj::ExactCountingAdapter>;

// Define error cases for weighted sketch types
CONTAIN_OVERLOAD_FAIL(wjRMFinal)
CONTAIN_OVERLOAD_FAIL(wjBFFinal)
#undef CONTAIN_OVERLOAD_FAIL

// Type for counting b-bit minhash with uint16_t counter (up to 65536)
using CBBMinHashType = mh::CountingBBitMinHasher<uint64_t, uint16_t>; // Is counting to 65536 enough for a transcriptome?

// Type for super minhash
using SuperMinHashType = mh::SuperMinHash<>;

// Enumeration of supported sketch types
enum Sketch: int {
    HLL,                // HyperLogLog
    BLOOM_FILTER,       // Bloom Filter
    RANGE_MINHASH,      // Range MinHash/KMV
    FULL_KHASH_SET,     // Full Hash Set
    COUNTING_RANGE_MINHASH, // Counting Range MinHash
    BB_MINHASH,         // B-bit MinHash
    BB_SUPERMINHASH,    // B-bit SuperMinHash
    COUNTING_BB_MINHASH,// Counting B-bit MinHash
    WIDE_HLL,           // Wide HyperLogLog
    HYPERMINHASH,       // HyperMinHash
    HMH = HYPERMINHASH  // Alias for HyperMinHash
};

// String names for each sketch type
static constexpr const char *const sketch_names [] {
    "HLL/HyperLogLog",
    "BF/BloomFilter", 
    "RMH/Range Min-Hash/KMV",
    "FHS/Full Hash Set",
    "CRHM/Counting Range Minhash",
    "BB/B-bit Minhash",
    "BBS/B-bit SuperMinHash",
    "CBB/Counting B-bit Minhash",
    "WHLL/Wide HLL",
    "HMH/HyperMinHash"
};

// Global arguments struct for sketch parameters
struct GlobalArgs {
    uint32_t weighted_jaccard_cmsize = 22;  // Count-min sketch size for weighted Jaccard
    uint32_t weighted_jaccard_nhashes = 10; // Number of hashes for weighted Jaccard
    uint32_t bbnbits = 16;                  // Number of bits for b-bit minhash
    uint32_t number_neighbors = 0;          // Number of neighbors (0 = disabled)
    size_t nperbatch = 16;                  // Number of sketches to process per batch
    bool exact_weighted = false;            // Use exact weighted calculations
    bool defer_hll_creation = false;        // Defer HLL creation
    void show() const {
        std::fprintf(stderr, "Global Arguments: %u wjcm, %u wjnh, %u bbits %u nn\n", weighted_jaccard_cmsize, weighted_jaccard_nhashes, bbnbits, number_neighbors);
    }
};

// Global variables
extern GlobalArgs gargs;
extern uint64_t global_hash_seed;

// Convert emission type to nearest neighbor type
INLINE static constexpr
NNType emt2nntype(EmissionType result_type) {
    switch(result_type) {
        case MASH_DIST: case FULL_MASH_DIST: case CONTAINMENT_DIST:
        case FULL_CONTAINMENT_DIST: case SYMMETRIC_CONTAINMENT_DIST:
            return DIST_MEASURE;
        case JI: case SIZES: case CONTAINMENT_INDEX: case SYMMETRIC_CONTAINMENT_INDEX:
        default:
            return SIMILARITY_MEASURE;
    }
    __builtin_unreachable();
    return SIMILARITY_MEASURE;
}
/**
 * @brief Converts EmissionType enum to string representation
 * @param result_type The EmissionType to convert
 * @return String representation of the EmissionType, or "ILLEGAL_EMISSION_FMT" if invalid
 */
static const char *emt2str(EmissionType result_type) {
    switch(result_type) {
        case MASH_DIST: return "MASH_DIST";
        case JI: return "JI"; 
        case SIZES: return "SIZES";
        case FULL_MASH_DIST: return "FULL_MASH_DIST";
        case FULL_CONTAINMENT_DIST: return "FULL_CONTAINMENT_DIST";
        case CONTAINMENT_INDEX: return "CONTAINMENT_INDEX";
        case CONTAINMENT_DIST: return "CONTAINMENT_DIST";
        case SYMMETRIC_CONTAINMENT_INDEX: return "SYMMETRIC_CONTAINMENT_INDEX";
        case SYMMETRIC_CONTAINMENT_DIST: return "SYMMETRIC_CONTAINMENT_DIST";
        default: break;
    }
    return "ILLEGAL_EMISSION_FMT";
}

template<typename SketchType>
struct FinalSketch {
    using final_type = SketchType;
};
// Macro to define template specialization for FinalSketch with single type parameter
// Extracts the final_type from the template parameter type
#define FINAL_OVERLOAD(x) \
template<> struct FinalSketch<x> { \
    using final_type = typename x::final_type;}

// Macro to define template specialization for FinalSketch with two type parameters
// Extracts the final_type from the template parameter types
#define FINAL_OVERLOAD2(x, y) \
template<> struct FinalSketch<x, y> { \
    using final_type = typename x, y::final_type;}
// Template specializations for various sketch types
FINAL_OVERLOAD(mh::CountingRangeMinHash<uint64_t>);     // Counting range min-hash
FINAL_OVERLOAD(mh::RangeMinHash<uint64_t>);             // Range min-hash
FINAL_OVERLOAD(BKHash64);                               // BK hash
FINAL_OVERLOAD(mh::BBitMinHasher<uint64_t>);            // b-bit min-hash
FINAL_OVERLOAD(WideHyperLogLogHasher<>);                // Wide HyperLogLog
FINAL_OVERLOAD(HyperLogLogHasher<>);                    // HyperLogLog
FINAL_OVERLOAD(SuperMinHashType);                       // Super min-hash
FINAL_OVERLOAD(CBBMinHashType);                         // CBB min-hash
FINAL_OVERLOAD(bf::bf_t);                               // Bloom filter
FINAL_OVERLOAD(hll::hll_t);                             // HyperLogLog
FINAL_OVERLOAD(khset64_t);                              // K-hash set
FINAL_OVERLOAD(HyperMinHash);                           // Hyper min-hash
FINAL_OVERLOAD(wj::WeightedSketcher<bf::bf_t>);         // Weighted Bloom filter
FINAL_OVERLOAD(wj::WeightedSketcher<hll::hll_t>);       // Weighted HyperLogLog
FINAL_OVERLOAD(wj::WeightedSketcher<khset64_t>);        // Weighted k-hash set
FINAL_OVERLOAD(wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>);  // Weighted counting range min-hash
FINAL_OVERLOAD(wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>);          // Weighted range min-hash
FINAL_OVERLOAD(wj::WeightedSketcher<BKHash64>);         // Weighted BK hash
FINAL_OVERLOAD(wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>);         // Weighted b-bit min-hash
FINAL_OVERLOAD(wj::WeightedSketcher<SuperMinHashType>); // Weighted super min-hash
FINAL_OVERLOAD(wj::WeightedSketcher<WideHyperLogLogHasher<>>);             // Weighted wide HyperLogLog
FINAL_OVERLOAD(wj::WeightedSketcher<HyperLogLogHasher<>>);                 // Weighted HyperLogLog
FINAL_OVERLOAD(wj::WeightedSketcher<CBBMinHashType>);   // Weighted CBB min-hash
FINAL_OVERLOAD(wj::WeightedSketcher<HyperMinHash>);     // Weighted hyper min-hash
FINAL_OVERLOAD2(wj::WeightedSketcher<bf::bf_t, wj::ExactCountingAdapter>); // Exact weighted Bloom filter
FINAL_OVERLOAD2(wj::WeightedSketcher<hll::hll_t, wj::ExactCountingAdapter>); // Exact weighted HyperLogLog
FINAL_OVERLOAD2(wj::WeightedSketcher<khset64_t, wj::ExactCountingAdapter>); // Exact weighted k-hash set
FINAL_OVERLOAD2(wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>, wj::ExactCountingAdapter>); // Exact weighted counting range min-hash
FINAL_OVERLOAD2(wj::WeightedSketcher<mh::RangeMinHash<uint64_t>, wj::ExactCountingAdapter>); // Exact weighted range min-hash
FINAL_OVERLOAD2(wj::WeightedSketcher<BKHash64, wj::ExactCountingAdapter>); // Exact weighted BK hash
FINAL_OVERLOAD2(wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>, wj::ExactCountingAdapter>); // Exact weighted b-bit min-hash
FINAL_OVERLOAD2(wj::WeightedSketcher<SuperMinHashType, wj::ExactCountingAdapter>); // Exact weighted super min-hash
FINAL_OVERLOAD2(wj::WeightedSketcher<WideHyperLogLogHasher<>, wj::ExactCountingAdapter>); // Exact weighted wide HyperLogLog
FINAL_OVERLOAD2(wj::WeightedSketcher<HyperLogLogHasher<>, wj::ExactCountingAdapter>); // Exact weighted HyperLogLog
FINAL_OVERLOAD2(wj::WeightedSketcher<CBBMinHashType, wj::ExactCountingAdapter>); // Exact weighted CBB min-hash
FINAL_OVERLOAD2(wj::WeightedSketcher<HyperMinHash, wj::ExactCountingAdapter>); // Exact weighted hyper min-hash
// Default template for sketch file suffixes - uses ".sketch" as default suffix
template<typename T>struct SketchFileSuffix {static constexpr const char *suffix = ".sketch";};

// Macro to define specializations for sketch file suffixes
// For each sketch type, defines 3 specializations:
// 1. Base type with given suffix
// 2. WeightedSketcher version with ".wj" + suffix
// 3. ExactCounting WeightedSketcher with ".wj.exact" + suffix
#define SSS(type, suf) \
    template<> struct SketchFileSuffix<type> {static constexpr const char *suffix = suf;};\
    template<> struct SketchFileSuffix<wj::WeightedSketcher<type>> {static constexpr const char *suffix = ".wj" suf;};\
    template<> struct SketchFileSuffix<wj::WeightedSketcher<type, wj::ExactCountingAdapter>> {static constexpr const char *suffix = ".wj.exact" suf;}

// Specializations for different sketch types with their corresponding file suffixes
SSS(mh::CountingRangeMinHash<uint64_t>, ".crmh");  // Counting range min-hash
SSS(mh::RangeMinHash<uint64_t>, ".rmh");           // Range min-hash
SSS(BKHash64, ".rmh");                             // BK hash (uses same suffix as RangeMinHash)
SSS(khset64_t, ".khs");                            // K-hash set
SSS(bf::bf_t, ".bf");                              // Bloom filter
SSS(mh::BBitMinHasher<uint64_t>, ".bmh");          // B-bit min-hash
SSS(WideHyperLogLogHasher<>, ".whll");             // Wide HyperLogLog
SSS(SuperMinHashType, ".bbs");                      // Super min-hash
SSS(CBBMinHashType, ".cbmh");                       // Counting BB min-hash
SSS(HyperMinHash, ".hmh");                         // Hyper min-hash
SSS(hll::hll_t, ".hll");                           // HyperLogLog
SSS(HyperLogLogHasher<>, ".hll");                  // HyperLogLog hasher (uses same suffix as hll_t)

// Clean up macros
#undef SSS
#undef FINAL_OVERLOAD

namespace detail {

// Structure to hold a path and its associated size
struct path_size {
    friend void swap(path_size&, path_size&); // Friend declaration for swap function
    std::string path;  // Path string
    size_t size;       // Size value

    // Constructor taking rvalue path and size
    path_size(std::string &&p, size_t sz): path(std::move(p)), size(sz) {}
    
    // Constructor taking const reference path and size 
    path_size(const std::string &p, size_t sz): path(p), size(sz) {}
    
    // Move constructor
    path_size(path_size &&o): path(std::move(o.path)), size(o.size) {}
    
    // Default constructor
    path_size(): size(0) {}
    
    // Move assignment operator
    path_size &operator=(path_size &&o) {
        std::swap(o.path, path);
        std::swap(o.size, size);
        return *this;
    }
};

// Swap function to exchange contents of two path_size objects
inline void swap(path_size &a, path_size &b) {
    std::swap(a.path, b.path);
    std::swap(a.size, b.size);
}

} // namespace detail

// Determines if a given emission type produces symmetric results
static constexpr bool is_symmetric(EmissionType result_type) {
    switch(result_type) {
        // These emission types are symmetric (A->B == B->A)
        case MASH_DIST: case JI: case SIZES:
        case FULL_MASH_DIST: case SYMMETRIC_CONTAINMENT_INDEX: case SYMMETRIC_CONTAINMENT_DIST:
            return true;

        // These emission types are asymmetric (A->B != B->A)
        case CONTAINMENT_INDEX: case CONTAINMENT_DIST: case FULL_CONTAINMENT_DIST:
        default: break;
    }
    return false;
}

// Type alias for HyperLogLogHasher
using HLLH = HyperLogLogHasher<>;

// Template struct to map sketch types to their corresponding enum values
template<typename T> struct SketchEnum;

// Specializations for various sketch types mapping to their enum values
template<> struct SketchEnum<hll::hll_t> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<HLLH> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<bf::bf_t> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<mh::RangeMinHash<uint64_t>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<BKHash64> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<mh::CountingRangeMinHash<uint64_t>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<mh::BBitMinHasher<uint64_t>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<CBBMinHashType> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<khset64_t> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<SuperMinHashType> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<WideHyperLogLogHasher<>> {static constexpr Sketch value = WIDE_HLL;};
template<> struct SketchEnum<HyperMinHash> {static constexpr Sketch value = HYPERMINHASH;};

// Specializations for weighted sketchers using default counting adapter
template<> struct SketchEnum<wj::WeightedSketcher<hll::hll_t>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HLLH>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<bf::bf_t>> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<BKHash64>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<CBBMinHashType>> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<khset64_t>> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<wj::WeightedSketcher<SuperMinHashType>> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<WideHyperLogLogHasher<>>> {static constexpr Sketch value = WIDE_HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HyperMinHash>> {static constexpr Sketch value = HYPERMINHASH;};

// Specializations for weighted sketchers using exact counting adapter
template<> struct SketchEnum<wj::WeightedSketcher<hll::hll_t, wj::ExactCountingAdapter>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HLLH, wj::ExactCountingAdapter>> {static constexpr Sketch value = HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<bf::bf_t, wj::ExactCountingAdapter>> {static constexpr Sketch value = BLOOM_FILTER;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::RangeMinHash<uint64_t>, wj::ExactCountingAdapter>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<BKHash64, wj::ExactCountingAdapter>> {static constexpr Sketch value = RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::CountingRangeMinHash<uint64_t>, wj::ExactCountingAdapter>> {static constexpr Sketch value = COUNTING_RANGE_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<mh::BBitMinHasher<uint64_t>, wj::ExactCountingAdapter>> {static constexpr Sketch value = BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<CBBMinHashType, wj::ExactCountingAdapter>> {static constexpr Sketch value = COUNTING_BB_MINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<khset64_t, wj::ExactCountingAdapter>> {static constexpr Sketch value = FULL_KHASH_SET;};
template<> struct SketchEnum<wj::WeightedSketcher<SuperMinHashType, wj::ExactCountingAdapter>> {static constexpr Sketch value = BB_SUPERMINHASH;};
template<> struct SketchEnum<wj::WeightedSketcher<WideHyperLogLogHasher<>, wj::ExactCountingAdapter>> {static constexpr Sketch value = WIDE_HLL;};
template<> struct SketchEnum<wj::WeightedSketcher<HyperMinHash, wj::ExactCountingAdapter>> {static constexpr Sketch value = HYPERMINHASH;};

/**
 * @brief Default implementation for setting estimation methods - no-op for non-HLL types
 * @param x The sketch object
 * @param estim Estimation method (unused)
 * @param jestim Joint estimation method (unused)
 */
template<typename T>
inline void set_estim_and_jestim(T &x, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {}

/**
 * @brief Specialization for HLL sketches to set estimation methods
 * @param h The HLL sketch
 * @param estim Estimation method to use
 * @param jestim Joint estimation method to use
 */
template<typename Hashstruct>
inline void set_estim_and_jestim(hll::hllbase_t<Hashstruct> &h, hll::EstimationMethod estim, hll::JointEstimationMethod jestim) {
    h.set_estim(estim);
    h.set_jestim(jestim);
}

// Forward declaration of construct function
template<typename T> inline T construct(size_t ssarg);

/**
 * @brief Template struct for constructing sketches based on whether they are weighted
 */
template<typename T, bool is_weighted> struct Constructor;

/**
 * @brief Constructs a sketch object of type T
 * @param ssarg Size argument for the sketch
 * @return Constructed sketch object
 */
template<typename T>
inline T construct(size_t ssarg) {
    Constructor<T, wj::is_weighted_sketch<T>::value> constructor;
    return constructor.create(ssarg);
}

/**
 * @brief Specialization for constructing non-weighted sketches
 */
template<typename T> struct Constructor<T, false> {
    static auto create(size_t ssarg) {
        return T(ssarg);
    }
};

/**
 * @brief Specialization for constructing BBitMinHasher sketches
 */
template<> struct Constructor<BBitMinHasher<uint64_t>, false> {
    static auto create(size_t ssarg) {
        return BBitMinHasher<uint64_t>(ssarg, gargs.bbnbits);
    }
};

/**
 * @brief Specialization for constructing weighted sketches
 */
template<typename T> struct Constructor<T, true> {
    static auto create(size_t ssarg) {
        using base_type = typename T::base_type;
        using cm_type = typename T::cm_type;
        return T(cm_type(16, gargs.weighted_jaccard_cmsize, gargs.weighted_jaccard_nhashes), construct<base_type>(ssarg));
    }
};
/**
 * @brief Specialization for constructing weighted BBitMinHasher sketches with an adapter
 * @tparam Adapter Template adapter type that wraps BBitMinHasher
 */
template<template<typename> typename Adapter> struct Constructor<Adapter<BBitMinHasher<uint64_t>>, true> {
    using Type = Adapter<BBitMinHasher<uint64_t>>;
    
    /**
     * @brief Creates a weighted BBitMinHasher sketch with the specified adapter
     * @param ssarg Size argument for the sketch
     * @return Constructed weighted sketch
     */
    static auto create(size_t ssarg) {
        using base_type = BBitMinHasher<uint64_t>;
        using cm_type = typename Adapter<BBitMinHasher<uint64_t>>::cm_type;
        return Type(cm_type(16, gargs.weighted_jaccard_cmsize, gargs.weighted_jaccard_nhashes), construct<base_type>(ssarg));
    }
};

/**
 * @brief Gets cardinality estimate for a sketch
 * @tparam T Type of sketch
 * @param x Sketch to get cardinality for
 * @return Estimated cardinality
 */
template<typename T>
inline double cardinality_estimate(T &x) {
    return x.cardinality_estimate();
}

// Specializations for different sketch types
template<> inline double cardinality_estimate(hll::hll_t &x) {return x.report();}
template<> inline double cardinality_estimate(mh::FinalBBitMinHash &x) {return x.est_cardinality_;}
template<> inline double cardinality_estimate(mh::FinalDivBBitMinHash &x) {return x.est_cardinality_;}
template<> inline double cardinality_estimate(sketch::HyperMinHash &x) {return x.getcard();}

/**
 * @brief Generates filename for a sketch based on parameters
 * @tparam SketchType Type of sketch
 * @param path Base path
 * @param sketch_p Sketch size parameter
 * @param wsz Window size
 * @param k Kmer size
 * @param csz Count size
 * @param spacing Spacing pattern
 * @param suffix Optional suffix
 * @param prefix Optional prefix
 * @param enct Encoding type
 * @return Generated filename
 */
template<typename SketchType>
/**
 * @brief Generates a filename for a sketch based on its parameters
 * @param path Base path for the file
 * @param sketch_p Sketch size parameter
 * @param wsz Window size
 * @param k Kmer size
 * @param csz Count size
 * @param spacing Spacing pattern string
 * @param suffix Optional suffix to add to filename
 * @param prefix Optional prefix path
 * @param enct Encoding type (BONSAI, NTHASH, or CYCLIC)
 * @return Generated filename string
 */
static inline std::string make_fname(const char *path, size_t sketch_p, int wsz, int k, int csz, const std::string &spacing,
                       const std::string &suffix="", const std::string &prefix="",
                       EncodingType enct=BONSAI) {
    // Start with prefix if provided
    std::string ret(prefix);
    if(ret.size()) ret += '/';
    
    // Extract filename from path
    {
        const char *p, *p2;
        p = (p = std::strchr(path, FNAME_SEP)) ? p + 1: path;
        if(ret.size() && (p2 = strrchr(p, '/'))) ret += std::string(p2 + 1);
        else                                     ret += p;
    }
    
    // Add window size
    ret += ".w";
    ret + std::to_string(std::max(csz, wsz));
    ret += ".";
    
    // Add kmer size and spacing
    ret += std::to_string(k);
    ret += ".spacing";
    ret += spacing;
    ret += '.';
    
    // Add encoding type
    ret += enct == BONSAI ? "": enct == NTHASH ? "nt.": "cyclic.";
    
    // Add optional suffix
    if(suffix.size()) {
        ret += "suf";
        ret += suffix;
        ret += '.';
    }
    
    // Add sketch size and type suffix
    ret += std::to_string(sketch_p);
    ret += SketchFileSuffix<SketchType>::suffix;
    return ret;
}

namespace us {
/**
 * @brief Base template for calculating union size between two sketches
 * @tparam T Type of sketch
 * @param a First sketch
 * @param b Second sketch 
 * @throws NotImplementedError if not specialized for type T
 */
template<typename T> inline double union_size(const T &a, const T &b) {
    throw NotImplementedError(std::string("union_size not available for type ") + __PRETTY_FUNCTION__);
}

/**
 * @brief Base template for calculating intersection size between two sketches
 * @tparam T Type of sketch
 * @param a First sketch
 * @param b Second sketch
 * @throws NotImplementedError if not specialized for type T
 */
template<typename T> inline double intersection_size(const T &a, const T &b) {
    throw NotImplementedError(std::string("intersection_size not available for type ") + __PRETTY_FUNCTION__);
}

/**
 * @brief Macro to declare specializations for union and intersection size calculations
 * @param type The sketch type to specialize for
 */
#define US_DEC(type) \
template<> inline double union_size<type> (const type &a, const type &b) { \
    return a.union_size(b); \
} \
template<> inline double intersection_size<type> (const type &a, const type &b) { \
    return a.intersection_size(b); \
}

// Declare specializations for RMFinal type
US_DEC(RMFinal)
// Declare specializations for HyperMinHash type
US_DEC(sketch::HyperMinHash)
// Declare specializations for CRMFinal type
US_DEC(CRMFinal)
// Declare specializations for khset64_t type
US_DEC(khset64_t)
#undef US_DEC

/**
 * @brief Specialization for calculating intersection size between two HLL sketches
 * @param h1 First HLL sketch
 * @param h2 Second HLL sketch
 * @return Intersection size estimate using inclusion-exclusion principle
 */
template<> inline double intersection_size<hll::hllbase_t<>>(const hll::hllbase_t<> &h1, const hll::hllbase_t<> &h2) {
    return std::max(0., h1.creport() + h2.creport() - h1.union_size(h2));
}

/**
 * @brief Specialization for calculating union size between two b-bit MinHash sketches
 * @param a First b-bit MinHash sketch
 * @param b Second b-bit MinHash sketch
 * @return Union size estimate using Jaccard index
 */
template<> inline double union_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
    return (a.est_cardinality_ + b.est_cardinality_ ) / (1. + a.jaccard_index(b));
}

/**
 * @brief Specialization for calculating intersection size between two b-bit MinHash sketches
 * @param a First b-bit MinHash sketch
 * @param b Second b-bit MinHash sketch
 * @return Intersection size estimate using Jaccard index
 */
template<> inline double intersection_size<mh::FinalBBitMinHash> (const mh::FinalBBitMinHash &a, const mh::FinalBBitMinHash &b) {
    const double ji = a.jaccard_index(b);
    return ji * (a.est_cardinality_ + b.est_cardinality_ ) / (1. + ji);
}
} // namespace us

/**
 * @brief Calculates symmetric containment between two sketches
 * @param x First sketch
 * @param y Second sketch 
 * @return Symmetric containment value = intersection size / min(size1, size2)
 */
template<typename T>
inline auto symmetric_containment_func(const T &x, const T &y) {
    auto tmp = set_triple(x, y);
    return tmp[2] / (std::min(tmp[0], tmp[1]) + tmp[2]);
}

/**
 * @brief Compares two sketches and returns a similarity/distance measure based on the specified emission type
 * @param lhs First sketch
 * @param rhs Second sketch
 * @param result_type Type of similarity/distance measure to compute
 * @param ksinv Inverse of k-mer size (1/k) for distance calculations
 * @return Float value representing the requested similarity/distance measure
 */
template<typename ST>
float result_cmp(const ST &lhs, const ST &rhs, EmissionType result_type, double ksinv) {
    double ret = std::numeric_limits<double>::infinity();
    switch(result_type) {
        case FULL_MASH_DIST: case MASH_DIST: case JI: {
            ret = similarity<const ST>(lhs, rhs);
            if(result_type == MASH_DIST) ret = dist_index(ret, ksinv);
            else if(result_type == FULL_MASH_DIST) ret = full_dist_index(ret, ksinv);
        } break;
        case SYMMETRIC_CONTAINMENT_DIST: case SYMMETRIC_CONTAINMENT_INDEX: case SIZES: case FULL_CONTAINMENT_DIST: case CONTAINMENT_INDEX: case CONTAINMENT_DIST: {
            const auto triple = set_triple(lhs, rhs);
            ret = triple[2];
            if(result_type == SYMMETRIC_CONTAINMENT_INDEX || result_type == SYMMETRIC_CONTAINMENT_DIST) {
                ret /= (std::min(triple[0], triple[1]) + triple[2]);
                if(result_type == SYMMETRIC_CONTAINMENT_DIST) ret = containment_dist(ret, ksinv);
            } else if(result_type == FULL_CONTAINMENT_DIST || result_type == CONTAINMENT_DIST || result_type == CONTAINMENT_INDEX) {
                ret /= (triple[0] + triple[1] + triple[2]);
                if(result_type == CONTAINMENT_DIST) ret = containment_dist(ret, ksinv);
                else if(result_type == FULL_CONTAINMENT_DIST) ret = full_containment_dist(ret, ksinv);
            } // else, result_type is (SIZES), and we return ret
        } break;
        default: __builtin_unreachable();
    }
    return static_cast<float>(ret);
}

/**
 * @brief Gets reference to HLL sketch from a sketch type
 * @param s Input sketch
 * @return Reference to underlying HLL sketch
 */
template<typename SketchType> inline hll::hll_t &get_hll(SketchType &s);

/**
 * @brief Template specialization for getting HLL reference from HLL type
 * @param s Input HLL sketch
 * @return Reference to input HLL sketch
 */
template<> inline hll::hll_t &get_hll<hll::hll_t>(hll::hll_t &s) {return s;}

/**
 * @brief Gets const reference to HLL sketch from a sketch type
 * @param s Input sketch
 * @return Const reference to underlying HLL sketch
 */
template<typename SketchType> inline const hll::hll_t &get_hll(const SketchType &s);

/**
 * @brief Template specialization for getting const HLL reference from HLL type
 * @param s Input HLL sketch
 * @return Const reference to input HLL sketch
 */
/**
 * @brief Template specialization for getting const HLL reference from HLL type
 * @param s Input HLL sketch
 * @return Const reference to input HLL sketch
 */
template<> inline const hll::hll_t &get_hll<hll::hll_t>(const hll::hll_t &s) {return s;}

/**
 * @brief Helper struct for parallel sketch estimation
 * @tparam SketchType Type of sketch being used
 */
template<typename SketchType>
struct est_helper {
    const Spacer                      &sp_;     // Spacer for k-mer generation
    const std::vector<std::string> &paths_;     // Input file paths
    std::mutex                         &m_;     // Mutex for thread safety
    const u64                          np_;     // Number of registers (2^np)
    const bool                      canon_;     // Whether to canonicalize k-mers
    void                            *data_;     // Optional auxiliary data
    std::vector<SketchType>         &hlls_;    // Vector of sketch objects
    kseq_t                            *ks_;    // FASTA/Q parser state
};

/**
 * @brief Helper function for parallel sketch estimation
 * @tparam SketchType Type of sketch being used
 * @tparam ScoreType Scoring scheme for k-mers
 * @param data_ Pointer to est_helper struct
 * @param index Index of current file
 * @param tid Thread ID
 */
template<typename SketchType, typename ScoreType=score::Lex>
void est_helper_fn(void *data_, long index, int tid) {
    est_helper<SketchType> &h(*(est_helper<SketchType> *)(data_));
    fill_lmers<ScoreType, SketchType>(h.hlls_[tid], h.paths_[index], h.sp_, h.canon_, h.data_, h.ks_ + tid);
}

/**
 * @brief Fill a sketch with k-mers from input files
 * @tparam SketchType Type of sketch being used
 * @tparam ScoreType Scoring scheme for k-mers
 * @param ret Output sketch to fill
 * @param paths Input file paths
 * @param k K-mer size
 * @param w Window size
 * @param spaces Spacing vector
 * @param canon Whether to canonicalize k-mers
 * @param data Optional auxiliary data
 * @param num_threads Number of threads to use
 * @param np Number of registers (2^np)
 * @param ks FASTA/Q parser state
 */
template<typename SketchType, typename ScoreType=score::Lex>
void fill_sketch(SketchType &ret, const std::vector<std::string> &paths,
              unsigned k, uint16_t w, const spvec_t &spaces, bool canon=true,
              void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=nullptr) {
    // Default to using all available threads if num_threads is negative.
    if(num_threads < 0) {
        num_threads = std::thread::hardware_concurrency();
        LOG_INFO("Number of threads was negative and has been adjusted to all available threads (%i).\n", num_threads);
    }
    const Spacer space(k, w, spaces);
    if(num_threads <= 1) {
        LOG_DEBUG("Starting serial\n");
        for(u64 i(0); i < paths.size(); fill_lmers<ScoreType, SketchType>(ret, paths[i++], space, canon, data, ks));
    } else {
        LOG_DEBUG("Starting parallel\n");
        std::mutex m;
        KSeqBufferHolder kseqs(num_threads);
        std::vector<SketchType> sketches;
        while(sketches.size() < (unsigned)num_threads) sketches.emplace_back(ret.clone());
        est_helper<SketchType> helper{space, paths, m, np, canon, data, sketches, kseqs.data()};
        kt_for(num_threads, &est_helper_fn<SketchType, ScoreType>, &helper, paths.size());
        auto &rhll = get_hll(ret);
        for(auto &sketch: sketches) rhll += get_hll(sketch);
    }

}

/**
 * @brief Create a new HLL sketch from input files
 * @tparam ScoreType Scoring scheme for k-mers
 * @param paths Input file paths
 * @param k K-mer size
 * @param w Window size
 * @param spaces Spacing vector
 * @param canon Whether to canonicalize k-mers
 * @param data Optional auxiliary data
 * @param num_threads Number of threads to use
 * @param np Number of registers (2^np)
 * @param ks FASTA/Q parser state
 * @param estim HLL estimation method
 * @param jestim Joint estimation method
 * @param clamp Whether to clamp estimates
 * @return New HLL sketch
 */
/**
 * @brief Create a new HLL sketch from input files
 * @tparam ScoreType Scoring scheme for k-mers (defaults to score::Lex)
 * @param paths Input file paths to process
 * @param k K-mer size
 * @param w Window size
 * @param spaces Spacing vector for spaced seeds
 * @param canon Whether to canonicalize k-mers (default: true)
 * @param data Optional auxiliary data (default: nullptr)
 * @param num_threads Number of threads to use (default: 1)
 * @param np Number of registers (2^np) (default: 23)
 * @param ks FASTA/Q parser state (default: false)
 * @param estim HLL estimation method (default: ERTL_MLE)
 * @param jestim Joint estimation method (default: ERTL_JOINT_MLE)
 * @param clamp Whether to clamp estimates (default: true)
 * @return New HLL sketch containing the k-mer set
 */
template<typename ScoreType=score::Lex>
hll::hll_t make_hll(const std::vector<std::string> &paths,
                unsigned k, uint16_t w, spvec_t spaces, bool canon=true,
                void *data=nullptr, int num_threads=1, u64 np=23, kseq_t *ks=false, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE, uint16_t jestim=hll::JointEstimationMethod::ERTL_JOINT_MLE, bool clamp=true) {
    // Create master HLL sketch with specified parameters
    hll::hll_t master(np, estim, (hll::JointEstimationMethod)jestim, 1, clamp);
    
    // Fill the sketch with k-mers from input paths
    fill_sketch<hll::hll_t, ScoreType>(master, paths, k, w, spaces, canon, data, num_threads, np, ks);
    
    return master;
}

/**
 * @brief Estimate cardinality of k-mer set from input files
 * @tparam ScoreType Scoring scheme for k-mers (defaults to score::Lex)
 * @param paths Input file paths to process
 * @param k K-mer size
 * @param w Window size
 * @param spaces Spacing vector for spaced seeds
 * @param canon Whether to canonicalize k-mers
 * @param data Optional auxiliary data (default: nullptr)
 * @param num_threads Number of threads to use (default: -1 for all available)
 * @param np Number of registers (2^np) (default: 23)
 * @param ks FASTA/Q parser state (default: nullptr)
 * @param estim HLL estimation method (default: ERTL_MLE)
 * @return Estimated cardinality of the k-mer set
 */
template<typename ScoreType=score::Lex>
u64 estimate_cardinality(const std::vector<std::string> &paths,
                            unsigned k, uint16_t w, spvec_t spaces, bool canon,
                            void *data=nullptr, int num_threads=-1, u64 np=23, kseq_t *ks=nullptr, hll::EstimationMethod estim=hll::EstimationMethod::ERTL_MLE) {
    // Create temporary HLL sketch and return its cardinality estimate
    auto tmp(make_hll<ScoreType>(paths, k, w, spaces, canon, data, num_threads, np, ks, estim));
    return tmp.report();
}

/**
 * @brief Compute pairwise distances between query and reference sketches
 * @tparam SketchType Type of sketch being used (e.g. HLL)
 * @param ofp Output file pointer to write results
 * @param hlls Array of sketches (references followed by queries)
 * @param inpaths Input file paths corresponding to sketches
 * @param use_scientific Whether to use scientific notation in output
 * @param k K-mer size used in sketches
 * @param result_type Type of distance/similarity measure to compute
 * @param emit_fmt Output format (BINARY, FULL_TSV, etc)
 * @param buffer_flush_size Buffer size for flushing output
 * @param nq Number of query sketches
 */
template<typename SketchType>
void partdist_loop(std::FILE *ofp, SketchType *hlls, const std::vector<std::string> &inpaths, const bool use_scientific, const unsigned k, const EmissionType result_type, EmissionFormat emit_fmt, const size_t buffer_flush_size,
                   size_t nq)
{
    // Calculate inverse of k for distance calculations
    const float ksinv = 1./ k;
    
    // Validate number of queries
    if(nq >= inpaths.size()) {
        UNRECOVERABLE_ERROR(ks::sprintf("Wrong number of query/references. (ip size: %zu, nq: %zu\n", inpaths.size(), nq).data());
    }
    
    // Calculate number of reference sketches
    size_t nr = inpaths.size() - nq;

#if TIMING
    auto start = std::chrono::high_resolution_clock::now();
#endif

    // Setup async writing and formatting futures
    std::future<void> write_future, fmt_future;
    
    // Initialize output buffers
    std::array<ks::string, 2> buffers;
    for(auto &b: buffers) b.resize(4 * nr);

    // Process each query sketch
    for(size_t qi = nr; qi < inpaths.size(); ++qi) {
        auto &hq = hlls[qi];
        
        // Allocate array for distances
        std::unique_ptr<float[]> arr(new float[nr]);
        
        // Calculate distances in parallel
        OMP_PFOR_DYN
        for(size_t j = 0; j < nr; ++j) {
            arr[j] = result_cmp(hlls[j], hq, result_type, ksinv);
        }

        // Handle different output formats
        switch(emit_fmt) {
            case BINARY:
                // Write binary output asynchronously
                if(write_future.valid()) write_future.get();
                write_future = std::async(std::launch::async, [arr=std::move(arr), nr,ofp]() {
                    if(unlikely(std::fwrite(arr.get(), sizeof(float), nr, ofp) != nr))
                    UNRECOVERABLE_ERROR("Error writing to binary file");
                });
                break;
                
            case UT_TSV: case UPPER_TRIANGULAR: default:
            case FULL_TSV:
                // Format and write TSV output asynchronously
                if(fmt_future.valid()) fmt_future.get();
                fmt_future = std::async(std::launch::async, [nr,qi,ofp,ind=qi-nr,&inpaths,use_scientific,arr=std::move(arr),&buffers,&write_future]() {
                    auto &buffer = buffers[qi & 1];
                    buffer += inpaths[qi];
                    for(size_t i = 0; i < nr; ++i)
                        buffer.sprintf("\t%g", arr[i]);
                    buffer.putc_('\n');
                    if(write_future.valid()) write_future.get();
                    write_future = std::async(std::launch::async, [ofp,&buffer]() {buffer.flush(::fileno(ofp));});
                });
                break;
        }
    }

    // Wait for pending operations to complete
    if(fmt_future.valid()) fmt_future.get();
    if(write_future.valid()) write_future.get();

#if TIMING
    auto end = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "partdist (%zu by %zu) took %gms\n", nr, nq, std::chrono::duration<double, std::milli>(end - start).count());
#endif
}
// Global executable name
static const char *executable = nullptr;

/**
 * @brief Get the executable name
 * @return String containing executable name or "unspecified" if not set
 */
static std::string get_executable() {
    return executable ? std::string(executable): "unspecified";
}

// Usage/help functions
void main_usage(char **argv);      // Display main program usage
void dist_usage(const char *arg);   // Display distance calculation usage
void sketch_usage(const char *arg); // Display sketching usage
void sketch_by_seq_usage(const char *arg); // Display sequence sketching usage
void flatten_usage();              // Display flatten command usage
void union_usage [[noreturn]] (char *ex); // Display union command usage

// Main command functions
int sketch_main(int argc, char *argv[]);        // Create sketches from input sequences
int fold_main(int argc, char *argv[]);          // Fold/compress sketches
int card_main(int argc, char *argv[]);          // Calculate cardinality estimates
int panel_main(int argc, char *argv[]);         // Process sketch panels
int dist_main(int argc, char *argv[]);          // Calculate distances between sketches
int print_binary_main(int argc, char *argv[]);  // Print binary format sketches
int mkdist_main(int argc, char *argv[]);        // Make distance matrices
int flatten_main(int argc, char *argv[]);       // Flatten sketch data structures
int hll_main(int argc, char *argv[]);           // HyperLogLog operations
int union_main(int argc, char *argv[]);         // Union multiple sketches
int view_main(int argc, char *argv[]);          // View sketch contents
int sketch_by_seq_main(int argc, char *argv[]); // Create sketches per sequence
int dist_by_seq_main(int argc, char *argv[]);   // Calculate distances per sequence
}

#endif /* DASHING_H__ */
