##############################################################################
##############################################################################

# Implementation of the Lloyd Algorithm for Clustering

import csv
import math
import random
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import f_oneway

with open('data.txt', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

for i in range(len(data)):
    for j in range(2):
        data[i][j] = float(data[i][j])


for element in data:
    element[0] = element[0]/10

def distance(point1, point2):
    line = 0
    for i in range(len(point1)):
        line += (point1[i] - point2[i])**2

    return math.sqrt(line)


def clusterConv(clusters):
    #clusters are stored in a dictionary, hence for each index it's corresponding cluster (index) is stored
    # here, this is converted into a structure where each cluster is stored with its corresponding data points
    values = []
    cluster = {}
    for x in clusters.values():
        values.append(x)

    clusterCopy = clusters.items()

    for i in range(len(values)):
        keys = [key  for (key, value) in clusterCopy if value == values[i]]
        cluster[values[i]] = list(keys)

    return(cluster)

def centerstoclusters(center, data):
    #after centers have been selected, assign each data point to the cluster corresponding to its nearest center
    clusters = {}
    
    for k in range(len(data)):
        distances = []
        for j in range(len(center)):
            distances.append(distance(center[j], data[k]))
        mindist = float('inf')
        for i in range(len(distances)):
            if distances[i] < mindist:
                mindist = distances[i]
                index = i
        clusters[k] = index
    
    Clusters = clusterConv(clusters)
            
    return Clusters

def centerofgravity(datapoint):
    center = 0
    for element in datapoint:
        center += element

    return center/len(datapoint)

def clusterstocenters(cluster):
    #after data points have been assigned to clusters, assign each cluster's center of gravity to be the clusters new center
    xCluster = []
    yCluster = []
    for element in cluster:
        xCluster.append(element[0])
        yCluster.append(element[1])
    
    x = centerofgravity(xCluster)
    y = centerofgravity(yCluster)
    
    center = [x, y]
    
    return center

def kMeansInitializer(k, data):
    #initialize k random centers
    
    centers = []

    while len(centers)  < k:
        anumber = random.randint(0, len(data))
        centers.append(data[anumber])

    return centers

def plotter(data, centers):

    df = pd.DataFrame(data, columns=['x', 'y'])
    
    sns.scatterplot(data = df, x=df["x"], y=df["y"])
    for element in centers:
        plt.scatter(x=element[0], y=element[1], color='r')
    plt.show()    
    
    return

def cluster_plotter(data, adict, centers):
    clusters = []
    for value in adict.values():
        values = []
        for index in value:
            values.append(data[index])
        clusters.append(values)
    for element in clusters:
        df = pd.DataFrame(element, columns=['x', 'y'])
        sns.scatterplot(data = df, x=df["x"], y=df["y"])
    for element in centers:
        plt.scatter(x=element[0], y=element[1], color='r')
    plt.show()
    return

def has_converged(old_center, new_center):
    return (sorted(old_center) == sorted(new_center))

def distancetocenter(center, data):
    distances = []
    for i in range(len(data)):
        distances.append(distance(center, data[i]))

    return distances
    
    
def SSE(centers, cluster, data):
    distances = []
    for i in range(len(centers)):
        indices = cluster[i]
        coordinates = []
        for j in range(len(indices)):
            coordinates.append(data[j])
            
        dist = distancetocenter(centers[i], coordinates)
        distances.append(dist)

    return f_oneway(*distances)[0]
        

def LloydAlgorithm(k, df):
    #initialize k random centers
    centers = kMeansInitializer(k, data)
    new_center = kMeansInitializer(k, data)
    
    while not has_converged(centers, new_center):
        
        centers = new_center

        
        CentToCl = centerstoclusters(centers, data)

        new_center = []
        for value in CentToCl.values():
            
            values = []
            for index in value:
                values.append(data[index])
            new_center.append(clusterstocenters(values))
            
    #plotter(data, new_center)
    cluster_plotter(data, CentToCl, new_center)
    
    #print(SSE(new_center, CentToCl, data))
        
    return 

#for i in range(2, 11):
#    print(LloydAlgorithm(i, data))

print(LloydAlgorithm(5, data))