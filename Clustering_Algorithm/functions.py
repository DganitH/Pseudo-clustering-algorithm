import re
import numpy as np


def load_file(file_name):
    """
    This function gets a file name and split the file into lines.
    """
    f = open(file_name, "r")
    text = f.read()
    return re.split('\n', text)


def replace_letter(string, i, let):
    """
    This function gets a string, a location- i and a letter- let.
    It returns the string after replacing the letter in location i with the letter let.
    """
    new_string = []
    for j in range(len(string)):
        if i != j:
            new_string += string[j]
        else:
            new_string += let
    list_to_str = ''.join([str(elem) for elem in new_string])
    return list_to_str


def insert_letter(string, i, let):
    """
    This function gets a string, a location- i and a letter- let.
    It returns the string after inserting the letter let in location i.
    """
    new_string = []
    for j in range(len(string)):
        if i != j:
            new_string += string[j]
        else:
            new_string += let
            new_string += string[j]
    listToStr = ''.join([str(elem) for elem in new_string])
    return listToStr


def delete_letter(string, i):
    """
    This function gets a string and a location- i.
    It returns the string after deleting the letter in location i.
    """
    new_string = []
    for j in range(len(string)):
        if i != j:
            new_string += string[j]
    listToStr = ''.join([str(elem) for elem in new_string[:7]])
    return listToStr


def create_sub_balls(indices):
    """
    This function gets a list of indices.
    For each index, it calculates its' substitution ball.
    The function returns a mapping from each index in the substitution balls to its' original index.
    """
    idx_dict = {}
    for idx in indices:
        for i in range(len(idx)):
            for let in ['A', 'C', 'G', 'T']:
                idx_dict[replace_letter(idx, i, let)] = idx
    return idx_dict


def create_ins_balls(indices):
    """
    This function gets a list of indices.
    For each index, it calculates its' insertion ball.
    The function returns a mapping from each index in the insertion balls to its' original index.
    """
    idx_dict = {}
    for idx in indices:
        for i in range(len(idx)):
            for let in ['A', 'C', 'G', 'T']:
                idx_dict[insert_letter(idx, i, let)] = idx
    return idx_dict


def create_del_balls(indices):
    """
    This function gets a list of indices.
    For each index, it calculates its' deletion ball.
    The function returns a mapping from each index in the deletion balls to its' original index.
    """
    idx_del_dict = {}
    for idx in indices:
        for i in range(len(idx)):
            idx_del_dict[delete_letter(idx, i)] = idx
    return idx_del_dict


def load_correct_division(file, index_len):
    """
    This function gets the file with the erroneous strands divided into groups according to their original strands,
    and the indices len.
    It returns two dictionaries:
    idx_to_strands - maps between the index of the original strand to the erroneous copies of this strand.
    strands_to_idx - maps between an erroneous strand to the index of the its' original strand.
    """
    text = load_file(file)

    idx_to_strands = {}
    strands_to_idx = {}

    prev_idx = text[0][:index_len]
    idx = prev_idx
    flag = 0

    for line in text:
        if len(line) and line != "*****************************":
            prev_idx = line[:index_len]
        if not len(line):
            flag = 0
        if flag == 1:
            idx_to_strands[idx] += [line]
            strands_to_idx[line] = idx
        if line == "*****************************":
            idx = prev_idx
            idx_to_strands[idx] = []
            flag = 1
    return idx_to_strands, strands_to_idx


def create_3mers_hist():
    """
    This function returns a hist of 3-mers, which means a dictionary in which the keys are sequences of 3 letters
    from {A,C,G,T}. All of the keys map to zero.
    """
    hist = {}
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            for k in ['A', 'C', 'G', 'T']:
                three_mer = i + j + k
                hist[three_mer] = 0
    return hist


def calc_3mers_hist(string):
    """
    This function gets a string and returns its' hist of 3-mers, which means a dictionary in which the keys are
    sequences of 3 letters from {A,C,G,T}. Each key maps to the amount of shows of this key in the string.
    """
    hist = create_3mers_hist()
    for i in range(len(string) - 2):
        three_mer = string[i] + string[i + 1] + string[i + 2]
        hist[three_mer] += 1
    return hist


def hist_to_array(hist):
    """
    This function gets a histogram of 3-mers and turns it into an array.
    """
    res_list = []
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            for k in ['A', 'C', 'G', 'T']:
                three_mer = i + j + k
                res_list += [hist[three_mer]]
    return np.array(res_list)


def best_idx(strand_hist, idx_to_hist):
    """
    This function gets a histogram of 3-mers of a strand and a dictionary that maps between an index and its' clusters'
    histogram of 3-mers. It returns the index of the closest histogram to the strands' histogram, according to norm L2.
    """
    best_idx = ""
    best_val = np.inf
    for idx in idx_to_hist:
        val = np.linalg.norm(strand_hist - idx_to_hist[idx])
        if val < best_val:
            best_val = val
            best_idx = idx
    return best_idx
