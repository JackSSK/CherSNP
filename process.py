class cov19_Processer(object):
    """docstring for Processer."""

    def __init__(self, patterns):
        super(cov19_Processer, self).__init__()
        self.patterns = []
        for ele in patterns:
            bounds = []
            for record in ele["patterns"]:
                temp = [record["start"], record["end"]]
                bounds.append(temp)
            bounds = sorted(bounds, key = lambda x: x[0],reverse=False)
            # Probabaly going need this later
            # genes = []
            # i = 0
            # next = 1
            # while True:
            #     cur = bounds[i]
            #     if next >= len(bounds):
            #         genes.append(cur)
            #         break
            #     else:suf = bounds[next]
            #     if cur[0] <= suf[0] <= cur[1] or suf[0] <= cur[0] <= suf[1]:
            #         if len(cur) >= len(suf):
            #             next += 1
            #         else:
            #             i = next
            #             next += 1
            #     elif cur[1] < suf[0]:
            #         genes.append(cur)
            #         i = next
            #         next = i + 1
            self.patterns.append({
                "id":ele["id"],
                "patterns":bounds
            })
