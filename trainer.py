import os
import re
from sklearn import svm
import gen_tools as t
import observer as obs

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

        # init part
        init_clf = self._fit("init")
        # t.encode_json([notation, obser])

        # term part
        term_clf = self._fit("term")

        # entCDS part
        entCDS_clf = self._fit("entCDS")

        # outCDS part
        outCDS_clf = self._fit("outCDS")


    # This function is to build classifier
    def _fit(self, name):
        clf = svm.SVC()
        correct = len(self.observ.subj[name]["correct"])
        wrong = len(self.observ.subj[name]["wrong"])
        notation = [1] * correct + [0] * wrong
        obser = self.observ.subj[name]["correct"] + self.observ.subj[name]["wrong"]
        clf.fit(obser,notation)
        # test here
        answer = clf.predict(obser)
        correct = 0
        for i in range(len(answer)):
            if answer[i] == notation[i]:
                correct += 1
        print(name, ': ', correct/len(answer))

        return clf
