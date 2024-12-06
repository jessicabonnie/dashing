#include "dashing.h"


namespace bns {

/**
 * @brief Flattens multiple distance matrices into a single binary file
 * 
 * @param fpaths Vector of file paths containing distance matrices to flatten
 * @param outpath Output path for the flattened binary file
 * @param k_values Vector of k-mer sizes used to generate the distance matrices
 * @return int 0 on success, 1 on allocation failure, 2 on file open failure
 *
 * This function takes multiple distance matrices generated with different k-mer sizes
 * and flattens them into a single binary file. The output format is:
 * - uint32_t: Number of k-mer sizes
 * - uint64_t: Number of matrix entries
 * - uint64_t: Number of sets/sequences
 * - unsigned[]: Array of k-mer sizes
 * - float[]: Flattened distance values
 */
int flatten_all(const std::vector<std::string> &fpaths, const std::string outpath, std::vector<unsigned> &k_values) {
    // Check for empty input
    if(fpaths.empty()) UNRECOVERABLE_ERROR("no fpaths, see usage.");
    
    const size_t nk = k_values.size();
    
    // Load all distance matrices
    std::vector<dm::DistanceMatrix<float>> dms;
    dms.reserve(nk);
    for(const auto &fp: fpaths)
        dms.emplace_back(fp.data());
        
    // Verify all matrices have same dimensions
    const uint64_t ne = dms.front().num_entries();
    assert(std::accumulate(dms.begin() + 1, dms.end(), true,
           [ne](bool val, const auto &x) {return val && x.num_entries() == ne;}));
           
    // Allocate memory for flattened output
    float *outp = static_cast<float *>(std::malloc(nk * ne * sizeof(float)));
    if(!outp) {
        std::fprintf(stderr, "Allocation of %zu bytes failed\n", size_t(nk * ne * sizeof(float))); 
        return 1;
    }

    // Process matrices in blocks for better cache utilization
    static constexpr uint64_t NB = 4096;
    #pragma omp parallel for
    for(size_t i = 0; i < ((NB - 1) + ne) / NB; ++i) {
        auto spos = i * NB, espos = std::min((i + 1) * NB, ne);
        auto destp = outp + spos * nk;
        // Copy values from each matrix into flattened array
        do {for(auto j = 0u; j < nk;*destp++ = dms[j++][spos]);} while(++spos < espos);
    }

    // Open output file
    std::FILE *ofp = fopen(outpath.data(), "wb");
    if(!ofp) return 2;
    
    uint64_t number_sets = fpaths.size();
    uint32_t numk = k_values.size();
#if 0
    gzread(ifp, &nk, sizeof(nk));
    std::vector<unsigned> k_values(nk);
    gzread(ifp, &number_entries, sizeof(number_entries));
    gzread(ifp, &number_sets, sizeof(number_sets));
    gzread(ifp, k_values.data(), k_values.size() * sizeof(unsigned));
#endif
    if(std::fwrite(&numk, sizeof(numk), 1, ofp) != 1) UNRECOVERABLE_ERROR("Failed to write nk");
    if(std::fwrite(&ne, sizeof(ne), 1, ofp) != 1) UNRECOVERABLE_ERROR("Failed to write num entries");
    std::fwrite(&number_sets, sizeof(number_sets), 1, ofp);
    std::fprintf(stderr, "Wrote %u nk and %zu num e and %zu num sets\n", numk, size_t(ne), size_t(number_sets));
    
    // Write k-mer sizes
    if(std::fwrite(k_values.data(), sizeof(unsigned), k_values.size(), ofp) != k_values.size()) 
        UNRECOVERABLE_ERROR("Wrong count of written values");
    
    // Write flattened distance values
    size_t nwritten;
    if((nwritten = std::fwrite(outp, sizeof(float), nk * ne, ofp)) != nk * ne) {
        UNRECOVERABLE_ERROR("Failed to write nk * ne");
    }
    std::fprintf(stderr, "Wrote %zu items (%zu * %zu) [%zu bytes]\n", size_t(nk * ne), nk, size_t(ne), size_t(nk * ne * 4));
    
    // Cleanup
    std::fclose(ofp);
    std::free(outp);
    return 0;
}

}
