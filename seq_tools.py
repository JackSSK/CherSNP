

def complementary(seq):
    ans = []
    for i in seq[::-1].upper():
        if i == 'A':
            ans.append('T')
        elif i == 'T':
            ans.append('A')
        elif i == 'G':
            ans.append('C')
        elif i == 'C':
            ans.append('G')
    return ''.join(ans)
