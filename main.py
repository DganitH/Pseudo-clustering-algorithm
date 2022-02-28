from functions import *
from results_functions import *

if __name__ == '__main__':

    """__________Preliminary work__________"""

    indices = load_file("450indices.txt")
    idx_len = len(indices[0])
    sub_idx_dict = create_sub_balls(indices)  # maps an index with at most one substitution error to the original index.
    del_idx_dict = create_del_balls(indices)  # maps an index with one deletion error to the original index.
    ins_idx_dict = create_ins_balls(indices)  # maps an index with one insertion error to the original index.

    # correct_idx_to_strands - maps between the index of the original strand to the erroneous copies of this strand.
    # correct_strands_to_idx - maps between an erroneous strand to the index of the its' original strand.
    correct_idx_to_strands, correct_strands_to_idx = load_correct_division("evyat.txt", idx_len)

    strands_with_errors = load_file("errors_shuffled.txt")

    """__________First division into clusters__________"""

    strand_to_idx = {}  # maps a strand to its' index.
    idx_to_strands = {}  # maps an index to its' strands.

    unmapped_strands = []  # a list of strands that weren't mapped in the first division.
    strands_to_candidates = {}  # maps a strand to some candidates for its' index.

    for idx in indices:
        idx_to_strands[idx] = []

    for strand in strands_with_errors:
        if strand[:idx_len] in indices:  # the index is in the original indices list.
            strand_to_idx[strand] = strand[:idx_len]
            idx_to_strands[strand[:idx_len]] += [strand]
        else:
            # the strands' index is in more than one dictionary.
            if int(strand[:idx_len + 1] in ins_idx_dict) + int(strand[:idx_len] in sub_idx_dict) + \
                    int(strand[:idx_len - 1] in del_idx_dict) >= 2:
                strands_to_candidates[strand] = []
                if strand[:idx_len + 1] in ins_idx_dict:
                    strands_to_candidates[strand] += [ins_idx_dict[strand[:idx_len + 1]]]
                if strand[:idx_len] in sub_idx_dict:
                    strands_to_candidates[strand] += [sub_idx_dict[strand[:idx_len]]]
                if strand[:idx_len - 1] in del_idx_dict:
                    strands_to_candidates[strand] += [del_idx_dict[strand[:idx_len - 1]]]
                continue
            else:
                # the strand's first idx_len letters belong to sub_idx_dict.
                if strand[:idx_len] in sub_idx_dict:
                    strand_to_idx[strand] = sub_idx_dict[strand[:idx_len]]
                    idx_to_strands[sub_idx_dict[strand[:idx_len]]] += [strand]
                else:
                    # the strand's first idx_len+1 letters belong to ins_idx_dict.
                    if strand[:idx_len + 1] in ins_idx_dict:
                        strand_to_idx[strand] = ins_idx_dict[strand[:idx_len + 1]]
                        idx_to_strands[ins_idx_dict[strand[:idx_len + 1]]] += [strand]
                    else:
                        # if we move this part from comment,
                        # we will be checking if the strand's first idx_len-1 letters belong to del_idx_dict.
                        """
                        if strand[:(idx_len - 1)] in del_idx_dict:
                            strand_to_idx[strand] = del_idx_dict[strand[:idx_len - 1]]
                            idx_to_strands[del_idx_dict[strand[:idx_len - 1]]] += [strand]
                        else:
                        """
                        unmapped_strands += [strand]

    "__________Testing the first division into clusters__________"
    # counting the strands that were divided correctly in the first division.
    correct = 0
    false = 0
    for strand in strand_to_idx:
        if strand_to_idx[strand] == correct_strands_to_idx[strand]:
            correct += 1
        else:
            false += 1

    # counting the strands that were divided correctly in each cluster.
    idx_to_correct_percentage = {}
    idx_correct = 0
    percentages = []
    for idx in indices:
        for strand in idx_to_strands[idx]:
            if strand_to_idx[strand] == correct_strands_to_idx[strand]:
                idx_correct += 1
        idx_to_correct_percentage[idx] = idx_correct / len(idx_to_strands[idx])
        percentages += [idx_to_correct_percentage[idx]]
        idx_correct = 0

    print_stat_first_level(len(strands_with_errors), len(strand_to_idx), correct, false, percentages)
    first_div_graph()

    """__________Second division into clusters__________"""

    idx_to_3mer = {}  # maps between an index and its' cluster 3mer histogram.

    # Creating each cluster its' cluster 3mer histogram.
    for idx in idx_to_strands:
        res = np.zeros((1, 64))
        for i in range(min(len(idx_to_strands[idx]), 10)):
            hist = calc_3mers_hist(idx_to_strands[idx][i])
            res += hist_to_array(hist)
        idx_to_3mer[idx] = res / min(len(idx_to_strands[idx]), 10)

    # Going over the strands that weren't mapped in the first level.

    # Strands that have some indices candidates. Their histograms will be compared only to the histograms of their
    # candidates indices' clusters.
    for strand in strands_to_candidates:
        strand_hist = hist_to_array(calc_3mers_hist(strand))
        partial_idx_to_3mer = {}
        for idx in strands_to_candidates[strand]:
            partial_idx_to_3mer[idx] = idx_to_3mer[idx]
        idx = best_idx(strand_hist, partial_idx_to_3mer)
        strand_to_idx[strand] = idx
        idx_to_strands[idx] += [strand]

    # All the of the rest strands. Their histograms will be compared to the histograms of all of the clusters.
    for strand in unmapped_strands:
        hist = hist_to_array(calc_3mers_hist(strand))
        idx = best_idx(hist, idx_to_3mer)
        strand_to_idx[strand] = idx
        idx_to_strands[idx] += [strand]

    "__________Testing the second division into clusters__________"

    # counting the strands that were divided correctly in the second division.
    correct_candidates = 0
    correct_no_candidates = 0
    for strand in strands_to_candidates:
        if strand_to_idx[strand] == correct_strands_to_idx[strand]:
            correct_candidates += 1
    for strand in unmapped_strands:
        if strand_to_idx[strand] == correct_strands_to_idx[strand]:
            correct_no_candidates += 1

    # counting the strands that were divided correctly in each cluster.
    idx_to_correct_percentage = {}
    percentages = []
    for idx in indices:
        idx_correct = 0
        for strand in idx_to_strands[idx]:
            if strand_to_idx[strand] == correct_strands_to_idx[strand]:
                idx_correct += 1
        idx_to_correct_percentage[idx] = idx_correct / len(idx_to_strands[idx])
        percentages += [idx_to_correct_percentage[idx]]

    print_stat_second_level(len(strands_with_errors), len(strands_to_candidates), len(unmapped_strands),
                            correct_candidates, correct_no_candidates, percentages)
    second_div_graph(percentages)
