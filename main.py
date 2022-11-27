from loading import load_directory
from kmers import *
from collections import Counter


def index_file(sequences, k):
    return Counter(list_kmers(sequences, k))


def list_file(sequences, k):
    lst = []

    # Construct the list of all the kmers of all the sequences
    for seq in sequences:
        lst.extend([min(x) for x in stream_kmers(seq, k)])

    # Sort the list
    lst.sort()

    return lst


def intersect_index(index, sequences, k):
    """ Create an index containing all the kmers of the input sequences
    """
    index_uniq = index.copy()
    query_uniq = 0
    intersection = 0

    for seq in sequences:
        for kmer, rkmer in stream_kmers(seq, k):
            minmer = min(kmer, rkmer)

            # Query not in index
            if minmer not in index_uniq:
                query_uniq += 1
            # Query in index => intersection
            else:
                intersection += 1
                index_uniq[minmer] -= 1
                if index_uniq[minmer] == 0:
                    del index_uniq[minmer]

    return sum(index_uniq.values()), intersection, query_uniq



def intersect_sorted_lists(lst1, lst2):
    idx1 = idx2 = 0
    # Dataset specific or intersection kmer counts
    A = inter = B = 0

    while (idx1 < len(lst1) and idx2 < len(lst2)):
        kmer1 = lst1[idx1]
        kmer2 = lst2[idx2]

        # Same kmer => intersection
        if kmer1 == kmer2:
            inter += 1
            idx1 += 1
            idx2 += 1
        # first list specific
        elif kmer1 < kmer2:
            A += 1
            idx1 += 1
        # second list specific
        else:
            B += 1
            idx2 += 1

    # Add remaining kmers of the non empty list
    A += len(lst1) - idx1
    B += len(lst2) - idx2

    return A, inter, B


def intersect_sorted_lists_2(lst1, lst2):
    idx1 = idx2 = 0
    # Dataset specific or intersection kmer counts
    A = inter = B = 0

    while (idx1 < len(lst1) and idx2 < len(lst2)):
        kmer1 = lst1[idx1]
        kmer2 = lst2[idx2]

        # Same kmer => intersection
        if kmer1 == kmer2:
            inter += 1
            idx1 += 1
            idx2 += 1
        # first list specific
        elif kmer1 < kmer2:
            A += 1
            idx1 += 1
        # second list specific
        else:
            B += 1
            idx2 += 1

    # Add remaining kmers of the non empty list
    A += len(lst1) - idx1
    B += len(lst2) - idx2

    return A, inter, B


def list_kmers(sequences, k):
    kmers = []
    for seq in sequences:
        kmers.extend([min(kmer, rkmer) for kmer, rkmer in stream_kmers(seq, k)])
    return kmers


def similarity(A, inter, B):
    # +1 added for pseudocount. Avoid divisions by 0
    A_similarity = inter / (inter + A + 1)
    B_similarity = inter / (inter + B + 1)

    return A_similarity, B_similarity


def jaccard(A, inter, B):
    return inter / (A + inter + B)

# Load all the files in a dictionary
files = load_directory("data")

k = 21

# Loading
indexes = {f:index_file(files[f], k) for f in files}
lists = {f:list_file(files[f], k) for f in files}
    
filenames = list(files.keys())

# is the list of hashed kmers
#print(lists[filenames[0]]) 
#print(sketch(100, stream_kmers, filenames[0], k))
print(intersect_sorted_lists(sketch(100, stream_kmers, filenames[0] , k), sketch(100, stream_kmers, filenames[1], k)))

if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")

    k = 21

    # Loading
    indexes = {f:index_file(files[f], k) for f in files}
    lists = {f:list_file(files[f], k) for f in files}
    
    filenames = list(files.keys())

    for i in range(len(files)):
        for j in range(i+1, len(files)):
            # Method 1 using index
            """print("index")
            print(indexes[filenames[i]], indexes[filenames[j]])
            A, inter, B = intersect_index(indexes[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))"""

            # Method 2 using sorted lists
            print("lists")
            #A, inter, B = intersect_sorted_lists(lists[filenames[i]], lists[filenames[j]])
            A, inter, B = intersect_sorted_lists(sketch(1000, stream_kmers, filenames[i] , k), sketch(1000, stream_kmers, filenames[j], k))
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
            print()
