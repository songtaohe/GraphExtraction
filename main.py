# 
# Please check https://github.com/mitroadmaps/roadtracer/blob/master/lib/discoverlib/mapextract.py 
# for the original graph extraction code. 
#

import scipy.ndimage
from scipy.ndimage.filters import gaussian_filter
import skimage.morphology
import os
from PIL import Image
import numpy
import numpy as np
from multiprocessing import Pool
import sys
from math import sqrt
import pickle
from postprocessing import graph_refine, connectDeadEnds
import cv2 
from douglasPeucker import simplifyGraph
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-input', action='store', dest='input', type=str,
                    help='Input segmentation file ', required =True)

parser.add_argument('-output', action='store', dest='output', type=str,
                    help='Output file prefix ', required =True)
                    
parser.add_argument('-threshold', action='store', dest='threshold', type=int,
                    help='Threshold for binarization (0-255) ', default = 127, required = False)

parser.add_argument('-gaussian_sigma', action='store', dest='gaussian_filter_sigma', type=str,
                    help='The sigma (pixels) of the gaussian filter. Set this to 0 to disable gaussian filter.', default = "3", required =False)                 

parser.add_argument('-grey_closing', action='store', dest='grey_closing', type=int,
                    help='The window size (pixels) for grey_closing', default = 6, required =False)                 

parser.add_argument('-connected_component_min_size', action='store', dest='connected_component_min_size', type=int,
                    help='Remove connected components whose total lengths are less than [connected_component_min_size] pixels', default = 50, required =False)                 

parser.add_argument('-spurs_min_size', action='store', dest='spurs_min_size', type=int,
                    help='Remove short edges (deadends) whose lengths are less than [spurs_min_size] pixels', default = 10, required =False)                 

parser.add_argument('-connect_deadends_threshold', action='store', dest='connect_deadends_threshold', type=int,
                    help='Deadends (terminal vertices) within [connect_deadends_threshold] pixels will be connected', default = 20, required =False)                 

parser.add_argument('-polyline_simplify_threshold', action='store', dest='polyline_simplify_threshold', type=int,
                    help='We use Douglas Peucker algorithm to simplify the polylines. A larger threshold gives more aggressive simplification', default = 5, required =False)                 


args = parser.parse_args()

def distance(a, b):
    return  sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)

def point_line_distance(point, start, end):
    if (start == end):
        return distance(point, start)
    else:
        n = abs(
            (end[0] - start[0]) * (start[1] - point[1]) - (start[0] - point[0]) * (end[1] - start[1])
        )
        d = sqrt(
            (end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2
        )
        return n / d

def rdp(points, epsilon):
    """
    Reduces a series of points to a simplified version that loses detail, but
    maintains the general shape of the series.
    """
    dmax = 0.0
    index = 0
    for i in range(1, len(points) - 1):
        d = point_line_distance(points[i], points[0], points[-1])
        if d > dmax:
            index = i
            dmax = d
    if dmax >= epsilon:
        results = rdp(points[:index+1], epsilon)[:-1] + rdp(points[index:], epsilon)
    else:
        results = [points[0], points[-1]]
    return results



# load image
im = scipy.ndimage.imread(args.input)
if len(im.shape) == 3:
	print 'warning: bad shape {}, using first channel only'.format(im.shape)
	im = im[:, :, 0]
im = numpy.swapaxes(im, 0, 1)

# apply gaussian filter 
sigma = float(args.gaussian_filter_sigma)
if sigma > 0:
    im = gaussian_filter(im, sigma=sigma)

# apply grey-scale closing 
windowSize = args.grey_closing
if windowSize > 1:
    im = scipy.ndimage.grey_closing(im, size=(6,6))


Image.fromarray(im).save(args.output+"_seg_before_binarization.png")

# binarization
im = im >= args.threshold

Image.fromarray(im).save(args.output+"_seg_after_binarization.png")

im = skimage.morphology.thin(im)
im = im.astype('uint8')

# extract a graph by placing vertices every THRESHOLD pixels, and at all intersections
vertices = []
edges = set()
def add_edge(src, dst):
	if (src, dst) in edges or (dst, src) in edges:
		return
	elif src == dst:
		return
	edges.add((src, dst))
point_to_neighbors = {}
q = []
while True:
	if len(q) > 0:
		lastid, i, j = q.pop()
		path = [vertices[lastid], (i, j)]
		if im[i, j] == 0:
			continue
		point_to_neighbors[(i, j)].remove(lastid)
		if len(point_to_neighbors[(i, j)]) == 0:
			del point_to_neighbors[(i, j)]
	else:
		w = numpy.where(im > 0)
		if len(w[0]) == 0:
			break
		i, j = w[0][0], w[1][0]
		lastid = len(vertices)
		vertices.append((i, j))
		path = [(i, j)]

	while True:
		im[i, j] = 0
		neighbors = []
		for oi in [-1, 0, 1]:
			for oj in [-1, 0, 1]:
				ni = i + oi
				nj = j + oj
				if ni >= 0 and ni < im.shape[0] and nj >= 0 and nj < im.shape[1] and im[ni, nj] > 0:
					neighbors.append((ni, nj))
		if len(neighbors) == 1 and (i, j) not in point_to_neighbors:
			ni, nj = neighbors[0]
			path.append((ni, nj))
			i, j = ni, nj
		else:
			if len(path) > 1:
				path = rdp(path, 2)
				if len(path) > 2:
					for point in path[1:-1]:
						curid = len(vertices)
						vertices.append(point)
						add_edge(lastid, curid)
						lastid = curid
				neighbor_count = len(neighbors) + len(point_to_neighbors.get((i, j), []))
				if neighbor_count == 0 or neighbor_count >= 2:
					curid = len(vertices)
					vertices.append(path[-1])
					add_edge(lastid, curid)
					lastid = curid
			for ni, nj in neighbors:
				if (ni, nj) not in point_to_neighbors:
					point_to_neighbors[(ni, nj)] = set()
				point_to_neighbors[(ni, nj)].add(lastid)
				q.append((lastid, ni, nj))
			for neighborid in point_to_neighbors.get((i, j), []):
				add_edge(neighborid, lastid)
			break

neighbors = {}

for edge in edges:

	nk1 = (vertices[edge[0]][1],vertices[edge[0]][0])
	nk2 = (vertices[edge[1]][1],vertices[edge[1]][0])
	
	if nk1 != nk2:
		if nk1 in neighbors:
			if nk2 in neighbors[nk1]:
				pass
			else:
				neighbors[nk1].append(nk2)
		else:
			neighbors[nk1] = [nk2]

		if  nk2 in neighbors:
			if nk1 in neighbors[nk2]:
				pass 
			else:
				neighbors[nk2].append(nk1)
		else:
			neighbors[nk2] = [nk1]

		
g = graph_refine(neighbors, isolated_thr = args.connected_component_min_size, spurs_thr=args.spurs_min_size)
g = connectDeadEnds(g, thr = args.connect_deadends_threshold)
g = simplifyGraph(g, e=args.polyline_simplify_threshold)


# save graph in pickle format 
pickle.dump(g, open(args.output+"graph.p", "wb"))

# save graph in txt format 
with open(args.output+"graph.txt", 'w') as f:
    nodeids = {}
    nid = 0
    for nloc, _ in g.items():
        nodeids[nloc] = nid 
        f.write('{} {}\n'.format(nloc[1], nloc[0]))
        nid+=1

    f.write('\n')

    for nloc, nei in g.items():
        nid1 = nodeids[nloc]
        for nn in nei:
            nid2 = nodeids[nn]
            f.write('{} {}\n'.format(nid1, nid2))


# visualization of the graph
dim = np.shape(im)
img = np.zeros((dim[0], dim[1]), dtype= np.uint8)
for nloc, nei in g.items():
    x1,y1 = int(nloc[1]), int(nloc[0])
    for nn in nei:
        x2,y2 = int(nn[1]), int(nn[0])
        cv2.line(img, (x1,y1), (x2,y2), (200), 3)

for nloc, nei in g.items():
    x1,y1 = int(nloc[1]), int(nloc[0])
    cv2.circle(img, (x1,y1), 5, (255), -1)

cv2.imwrite(args.output+"graph_vis.png", img)















