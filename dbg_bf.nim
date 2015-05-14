# not explained well in docs:
# 1 - how to return array of strings (e.g., get k-mers)
# 2 - access global table/array in function (used in get_rc)

import tables
import bloom

const
    bases = ['A', 'C', 'G', 'T']
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}.toTable()
    k = 27
    fp = 0.01
    j = 10
    read_len = 100

proc get_kmers(r: string, k: int, kmers: var openarray[string] ) =
    var i = 0
    while i < len(r)-k+1:
        kmers[i] = r[i .. i+k-1]
        i = i+1

proc reverse(s: var string): string =
    result = newString(len(s))
    for i in 0 .. s.high:
        result[i] = s[s.high - i]

proc get_rc(dna: var string): string = 
    result = reverse(dna) 
    for i in 0 .. dna.high:
        result[i] = complements[result[i]]
   

when isMainModule:
    var read = "ACGTTCGTTTGACACTTCGTTTGTCGTTTGGTTCGTTGTTCGTT"
    var kmers: array[0..read_len-k+1, string]
    # get_kmers(read, k, kmers)
    # for i,value in @kmers:
    #     if value!=nil:
    #         echo i, " ", value


    echo reverse(read)
    echo get_rc(read)
    echo read # original isn't changed
    
    var bf = initialize_bloom_filter(capacity = 10_000, error_rate = fp)

    # iterate over lines
    # if line starts with ">", skip
    var f_hand = open("chr20.c10.reads.head")
    for line in f_hand.lines:
        get_kmers(line,k,kmers)
        for i,value in @kmers:
            if value!=nil:
                echo i, " ", value

        # continue
        # process(line)
    i.close()

# # if ACGT line get k-mers as list
# # store first, last, last - 10 k-mers 
# # put all in BF
# # import bloom
# echo(bf)                                               # Get characteristics of the Bloom filter
# echo(bf.lookup("An element not in the Bloom filter"))  # Prints 'false'
# bf.insert("Here we go...")
# assert(bf.lookup("Here we go..."))


# def get_kmers(r,k):
#     return [r[i:i+k] for i in range(len(r)-k+1)]

# def get_rc(dna):
#     rev = reversed(dna)
#     return "".join([complements[i] for i in rev])
    
# def get_canons(kmers):
#     return [min(x,get_rc(x)) for x in kmers]