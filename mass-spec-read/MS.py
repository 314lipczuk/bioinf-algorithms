from collections import defaultdict
mz_values = [
    129.102, 175.118, 226.117, 244.129, 338.182, 355.160, 373.170, 
    460.208, 475.242, 559.779, 574.308, 624.329, 681.342, 702.367, 
    801.438, 902.478, 989.511, 1118.556
]

amino_acid_dict = {
    71.03711: "A",
    156.10111: "R",
    114.04293: "N",
    115.02694: "D",
    103.00919: "C",
    129.04259: "E",
    128.05858: "Q",
    57.02146: "G",
    137.05891: "H",
    113.08406: "I",
    113.08406: "L",
    128.09496: "K",
    131.04049: "M",
    147.06841: "F",
    97.05276: "P",
    87.03203: "S",
    101.04768: "T",
    186.07931: "W",
    163.06333: "Y",
    99.06841: "V"
}

vals = amino_acid_dict.keys()

#                           from        to
distance_range_to_check = [min(vals), max(vals)]

def possible_matches(a):
    result = set()
    for aa in vals:
        if abs(aa-a) < 0.02:
            result.add(amino_acid_dict[aa])
    return result

result = defaultdict(set)

d_n = []
for i in range(len(mz_values)-1) :
    peak = mz_values[i]
    candidates = [v for v in mz_values if peak < v and abs(peak - v) >= distance_range_to_check[0] and abs(peak - v) <= distance_range_to_check[1]]
    for c in candidates:
        d = abs(c - peak)
        matches = possible_matches(d)
        result[(peak, d)] = result[(peak, d)].union(matches)

for pair in result.keys():
    if len(result[pair]) > 0:
        print(f"{pair} -> {str(result[pair])}")
            
            
    

    
