from loading import load_directory
from kmers import stream_kmers, kmer2str



def similarity(A, inter, B):
    return inter / (A+inter), inter / (B+inter)


def jaccard(A, inter, B):
    return inter / (A + B +inter)


def intersection(A, B) :
    dict_1, dict_2 = {}, {}
    inter_kmers = []
    inter_A, inter_B = 0, 0
    inter_kmers = set()
    
    for kmer in A :
        if kmer not in dict_1 :
            dict_1[kmer] = 1
        else :
            dict_1[kmer] += 1
    
    for kmer in B :
        if kmer in dict_1 :
            inter_kmers.add(kmer)
            if kmer not in dict_2 :
                dict_2[kmer] = 1
            else :
                dict_2[kmer] += 1
                
    inter_kmers = set(inter_kmers)
    
    for kmer in inter_kmers :
        inter_A += dict_1[kmer]
        inter_B += dict_2[kmer]
    
    return inter_A, inter_B


def my_method(A, B) :
    inter_A, inter_B = intersection(A, B)
    print(inter_A, inter_B)
    return len(A)-inter_A, inter_A+inter_B, len(B)-inter_B

if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    filenames = list(files.keys())
    
    # Stream
    streamed_kmers = []
    for filename in filenames :
        kmers = []
        for i in files[filename] :
            kmers.extend(stream_kmers(i, k))
        files[filename] = kmers
        
    
    for i in range(len(files)):
        for j in range(i+1, len(files)):

            A, inter, B = my_method(files[filenames[i]], files[filenames[j]])
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
