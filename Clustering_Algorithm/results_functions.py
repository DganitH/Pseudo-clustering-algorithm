import statistics
import matplotlib.pyplot as plt
import numpy as np


def print_stat_first_level(strands_num, divided_strands_num, correct_num, false_num, percentages):
    """
    This function is showing some statistics after the first division.
    """
    print("__________Testing the first division into clusters__________")
    print("Total number of strands:", strands_num)
    print("The part of strands that were divided in the first division:",
          divided_strands_num / strands_num)
    print("The part of strands that were divided correctly in the first division (from the divided strands):",
          correct_num / divided_strands_num)
    print("The part of strands that were divided falsely in the first division (from the divided strands):",
          false_num / divided_strands_num)

    print("The minimum part of correct strands in a cluster:", min(percentages))
    print("The maximum part of correct strands in a cluster:", max(percentages))
    print("The mean of the part of correct strands in a cluster:", statistics.mean(percentages))
    print("The standard deviation of the part of correct strands in a cluster:", statistics.pstdev(percentages))
    print("The median of the part of correct strands in a cluster:", statistics.median(percentages))


def print_stat_second_level(strands_num, strands_to_candidates_num, unmapped_strands_num, correct_candidates_num,
                            correct_no_candidates_num, percentages):
    """
    This function is showing some statistics after the second division.
    """
    print("__________Testing the second division into clusters__________")
    print("The part of strands that were divided in the second division:",
          (strands_to_candidates_num + unmapped_strands_num) / strands_num)
    print("The amount of strands with candidates:", strands_to_candidates_num)
    print("The amount of strands with no candidates:", unmapped_strands_num)
    print("The part of strands that were divided correctly in the second division (from the divided strands):",
          (correct_candidates_num + correct_no_candidates_num) / (strands_to_candidates_num + unmapped_strands_num))
    print("The part of strands that were divided correctly from the strands with candidates:",
          correct_candidates_num / strands_to_candidates_num)
    print("The part of strands that were divided correctly from the strands with no candidates:",
          correct_no_candidates_num / unmapped_strands_num)

    print("The minimum part of correct strands in a cluster:", min(percentages))
    print("The maximum part of correct strands in a cluster:", max(percentages))
    print("The mean of the part of correct strands in a cluster:", statistics.mean(percentages))
    print("The standard deviation of the part of correct strands in a cluster:", statistics.pstdev(percentages))
    print("The median of the part of correct strands in a cluster:", statistics.median(percentages))


def first_div_graph():
    """
    This function is showing the graph of the part of correct strands in a cluster after the first division.
    """
    perc1 = np.load('perc1.npy')
    perc2 = np.load('perc2.npy')

    bins = np.linspace(0, 1, 100)
    plt.hist(perc1, bins, alpha=0.5, label='including deletions dictionary')
    plt.hist(perc2, bins, alpha=0.5, label='not including deletions dictionary')
    plt.legend(loc='upper left')
    plt.grid()
    plt.title('The part of correct strands in a cluster after the first division.')
    plt.ylabel('amount of clusters')
    plt.xlabel('part of correct strands in a cluster')
    plt.show()


def second_div_graph(percentages):
    """
    This function is showing the graph of the part of correct strands in a cluster after the second division.
    """
    bins = np.linspace(0, 1, 100)
    plt.hist(percentages, bins, alpha=0.5)
    plt.grid()
    plt.title('The part of correct strands in a cluster after the second division.')
    plt.ylabel('amount of clusters')
    plt.xlabel('part of correct strands in a cluster')
    plt.show()
