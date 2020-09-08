import os
import re
from sklearn import svm
# from itertools import chain
import gen_tools as t
import observer as obs

# import feature as fea
# Nothing here
# X = [[1,1,1], [2,2,2]]
# y = [0,1]
# clf = svm.SVC()
# clf.fit(X,y)
# print(clf.predict([[2,2,1]]))

class Trainer:
    def __init__(self, seq_file, gff_file):
        self.observ = obs.Observer(seq_file, gff_file)
        # test 
        t.encode_json([self.observ.dict, self.observ.subj])
