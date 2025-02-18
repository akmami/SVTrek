#include "refinement.h"

// /**
//  * consensus() picks the most common (or "peak") value in a range of integers
//  */
// int consensus(std::vector<int> &lengths, params &_params) {
//     int consensus_val = -1;
//     int max_count = _params.consensus_min_count - 1;

//     std::sort(lengths.begin(), lengths.end());
//     for (size_t i = 0; i < lengths.size(); i++) {
//         int count = 1;
//         // Count how many values fall within consensus_interval of lengths[i]
//         for (size_t j = i + 1; j < lengths.size() && lengths[j] <= lengths[i] + _params.consensus_interval; j++) {
//             count++;
//         }
//         if (count > max_count) {
//             max_count = count;
//             consensus_val = lengths[i];
//         }
//     }
//     return consensus_val;
// }

// /**
//  * consensus_pos() is similar, but picks the best cluster *and*
//  * tries to be closest to imprecise_pos
//  */
// int consensus_pos(std::vector<int> &locations, int imprecise_pos, params &_params) {
//     int consensus_val = -1;
//     int max_count = _params.consensus_min_count - 1;
//     int distance = INT_MAX;

//     std::sort(locations.begin(), locations.end());
//     for (size_t i = 0; i < locations.size(); i++) {
//         int count = 1;
//         for (size_t j = i + 1; j < locations.size() && locations[j] <= locations[i] + _params.consensus_interval; j++) {
//             count++;
//         }
//         int candidate = locations[i];
//         bool better_count = (count > max_count);
//         bool same_count_but_closer = (count == max_count && std::abs(imprecise_pos - candidate) < distance);

//         if (better_count || same_count_but_closer) {
//             max_count = count;
//             consensus_val = candidate;
//             distance = std::abs(imprecise_pos - candidate);
//         }
//     }
//     return consensus_val;
// }

// /**
//  * refine_start()
//  */
// void refine_start(int chrom, int outer_start, int inner_start, int imprecise_pos, result &res, params &_params, thread_data &_thread_data, svtype type) {
//     // if we don't use 'type', silence the warning
//     (void)type;

//     std::vector<int> start_positions;
//     bam1_t *aln = bam_init1();

//     // Query the region from [outer_start - wider_interval] to [inner_start + narrow_interval]
//     hts_itr_t *iter = sam_itr_queryi(_thread_data.htslib.bam_file_index,
//                                      chrom - 1,
//                                      outer_start - _params.wider_interval,
//                                      inner_start + _params.narrow_interval);

//     if (iter) {
//         while (sam_itr_next(_thread_data.htslib.fp_in, iter, aln) > 0) {
//             int reference_pos = aln->core.pos;
//             uint32_t *cigar   = bam_get_cigar(aln);

//             // n_cigar is unsigned, so use uint32_t in the loop
//             for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
//                 if (bam_cigar_op(cigar[i]) == __CIGAR_DELETION &&
//                     reference_pos >= outer_start &&
//                     reference_pos <= inner_start)
//                 {
//                     // +1 because positions are 1-based
//                     start_positions.push_back(reference_pos + 1);
//                 }
//                 // Advance reference_pos if the CIGAR operation consumes reference
//                 if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION &&
//                     bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP)
//                 {
//                     reference_pos += bam_cigar_oplen(cigar[i]);
//                 }
//                 if (reference_pos > inner_start)
//                     break;
//             }
//         }
//         sam_itr_destroy(iter);
//     }
//     bam_destroy1(aln);

//     // pick the best refined start
//     res.refined_start = consensus_pos(start_positions, imprecise_pos, _params);
// }

// /**
//  * refine_end()
//  */
// void refine_end(int chrom, int inner_end, int outer_end, int imprecise_pos, result &res, params &_params, thread_data &_thread_data, svtype type) {
//     (void)type;

//     std::vector<int> end_positions;
//     bam1_t *aln = bam_init1();

//     hts_itr_t *iter = sam_itr_queryi(_thread_data.htslib.bam_file_index,
//                                      chrom - 1,
//                                      inner_end - _params.wider_interval,
//                                      outer_end + _params.narrow_interval);

//     if (iter) {
//         while (sam_itr_next(_thread_data.htslib.fp_in, iter, aln) > 0) {
//             int reference_pos = aln->core.pos;
//             uint32_t *cigar   = bam_get_cigar(aln);

//             for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
//                 if (bam_cigar_op(cigar[i]) == __CIGAR_DELETION &&
//                     reference_pos >= inner_end &&
//                     reference_pos <= outer_end)
//                 {
//                     end_positions.push_back(reference_pos + 1);
//                 }
//                 if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION &&
//                     bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP)
//                 {
//                     reference_pos += bam_cigar_oplen(cigar[i]);
//                 }
//                 if (reference_pos > outer_end)
//                     break;
//             }
//         }
//         sam_itr_destroy(iter);
//     }
//     bam_destroy1(aln);

//     // pick the best refined end
//     res.refined_end = consensus_pos(end_positions, imprecise_pos, _params);
// }

// /**
//  * refine_point() - used for insertions or single-point refinement
//  */
// void refine_point(int chrom, int start, int end, int imprecise_pos,
//                   result &res,
//                   params &_params,
//                   thread_data &_thread_data,
//                   svtype type)
// {
//     (void)type;

//     std::vector<int> positions;
//     bam1_t *aln = bam_init1();

//     hts_itr_t *iter = sam_itr_queryi(_thread_data.htslib.bam_file_index,
//                                      chrom - 1,
//                                      start - _params.wider_interval,
//                                      end + _params.narrow_interval);

//     if (iter) {
//         while (sam_itr_next(_thread_data.htslib.fp_in, iter, aln) > 0) {
//             int reference_pos = aln->core.pos;
//             uint32_t *cigar   = bam_get_cigar(aln);

//             for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
//                 if (bam_cigar_op(cigar[i]) == __CIGAR_INSERTION &&
//                     reference_pos >= start &&
//                     reference_pos <= end)
//                 {
//                     positions.push_back(reference_pos + 1);
//                 }
//                 if (bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
//                     reference_pos += bam_cigar_oplen(cigar[i]);
//                 }
//                 if (reference_pos > end)
//                     break;
//             }
//         }
//         sam_itr_destroy(iter);
//     }
//     bam_destroy1(aln);

//     res.refined_start = consensus_pos(positions, imprecise_pos, _params);
// }


// result deletion(int chrom,
//                 int outer_start,
//                 int inner_start,
//                 int inner_end,
//                 int outer_end,
//                 int imprecise_pos,
//                 int imprecise_end,
//                 params &_params,
//                 thread_data &_thread_data)
// {

//     result res;
//     if ((imprecise_end - imprecise_pos) < __CI_MIN_FILTER_LENGTH ||
//         (imprecise_end - imprecise_pos) * _params.ci_max_length >
//         (inner_start - outer_start))
//     {
//         refine_start(chrom, outer_start, inner_start, imprecise_pos,
//                      res, _params, _thread_data, DELETION);
//     }
//     if ((imprecise_end - imprecise_pos) < __CI_MIN_FILTER_LENGTH ||
//         (imprecise_end - imprecise_pos) * _params.ci_max_length >
//         (outer_end - inner_end))
//     {
//         refine_end(chrom, inner_end, outer_end, imprecise_end,
//                    res, _params, _thread_data, DELETION);
//     }
//     return res;
// }

// result insertion(int chrom,
//                  int outer_start,
//                  int inner_start,
//                  int imprecise_pos,
//                  params &_params,
//                  thread_data &_thread_data)
// {
//     result res;
//     refine_point(chrom, outer_start, inner_start, imprecise_pos,
//                  res, _params, _thread_data, INSERTION);
//     return res;
// }

// result inversion(int chrom,
//                  int outer_start,
//                  int inner_start,
//                  int inner_end,
//                  int outer_end,
//                  int imprecise_pos,
//                  int imprecise_end,
//                  params &_params,
//                  thread_data &_thread_data)
// {
//     result res;
//     // Inversion might refine both ends similarly
//     if ((imprecise_end - imprecise_pos) < __CI_MIN_FILTER_LENGTH ||
//         (imprecise_end - imprecise_pos) * _params.ci_max_length >
//         (inner_start - outer_start))
//     {
//         refine_point(chrom, outer_start, inner_start, imprecise_pos,
//                      res, _params, _thread_data, INSERTION);
//     }
//     if ((imprecise_end - imprecise_pos) < __CI_MIN_FILTER_LENGTH ||
//         (imprecise_end - imprecise_pos) * _params.ci_max_length >
//         (outer_end - inner_end))
//     {
//         refine_point(chrom, inner_end, outer_end, imprecise_end,
//                      res, _params, _thread_data, INSERTION);
//     }
//     return res;
// }
