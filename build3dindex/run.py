"""
 * date: 2024-04-26
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  Obtain the 3d indexes for all datasets.
 """

import os
testObjFile = '/home/zoe/Documents/Experiment-viewshed/Q-view/build3dindex/test_dataset/paper_area500/'#absolute path

for root, dirs, files in os.walk(testObjFile):
    #print("Current directory:", root)
    for testObjFile in files:
        #print("Found file:", testObjFile)
        path = os.path.join(root, testObjFile)
        if (len(path.split('.'))>1) and (path.split('.')[1]=="obj"):
            #print("Processing .obj file:", path)
            p=path.split('/')
            cmd="./Build3dindex --input "+ path + " --output "+"./test_indexes/"+p[len(p)-2]+".3idx"
            print(cmd)
            os.system(cmd)
            

