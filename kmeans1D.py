import numpy as np


def KCenter1D(y,K):

	# Description
	# one-dimensional cluster algorithm implemented in Python
	# y is input one-dimensional vector and K stands for the cluster level
	# all vectors in this program is considered starting at position 1, position 0 is not used.
	# Input:
	# y -- a vector of numbers, not necessarily sorted
	# K -- the number of clusters expected

	# Pre-Conditions:
	# K <= |set(y)|
	# |set(y)| > 1
	if K>len(set(y)) or len(set(y))==1:
		return {'Cluster':[1]*len(y),'Centers':y,'Size':[1]*K}

	x = sorted(y)
	x = [None]+x   # 0 index is not used

	N = len(x) - 1 # N: is the size of input vector, -1 for None item

	#D - Matrix rows 1 to K, columns 1 to N
	#B - Matrix rows 1 to K, columns 1 to N
	D = np.array([[None]*(N+1)]*(K+1) ,dtype=np.float64)
	B = np.array([[None]*(N+1)]*(K+1) ,dtype=np.float64)

	#Set for all K Column i:noOfItem = 1 to 0
	for p in range(1,K+1):
			D[p][1] = 0.0
			B[p][1] = 1


	mean_x1     = 0.0     #mean of rightmost-cluster when k=1
	mean_xj		 = 0.0     #mean of rightmost-cluster when k>1

	d 					= 0.0     #SOS from mean for the rightmost-cluster

	for k in range(1,K+1):

		mean_x1 = x[1]

		for i in range(2,N+1):

			# when k = 1, D[1][i] = d(1,i)
			if (k == 1):

				D[1][i] = D[1][i-1] + (float(i-1)/float(i))*(x[i]-mean_x1)**2
				mean_x1 = (float64(i-1)*mean_x1 + x[i])/float64(i)
				B[1][i] = 1

			#when k > 1, D[k][i]  = min ( D[k-1][j-1] + d(j,i) ) for j>=1
			else:


				d = 0   			#Initiallly no item in rightmost-cluster
				mean_xj = 0	 #Initiallly no item in rightmost-cluster
				D[k][i] = -1  #value of D[k][i] is not set

				for j in reversed(range(1,i+1)):

					#One-by-one add item (j,i) to rightmost-cluster
					d += (float(i - j) / float(i-j+1)) * (x[j] - mean_xj) **2   #d(j,i)
					mean_xj = (x[j] + float(i-j)*mean_xj) / float(i-j+1)				#mean(j,i)

					#initialize D[k][i] with D[k][i] = D[k-1][i-1] + d
					if (D[k][i] == -1):
						if (j == 1):
							D[k][i] = d
							B[k][i] = float(j)
						else:
							D[k][i] = d + D[k-1][j-1]
							B[k][i] = float64(j)
					#If already initialized find if D[k-1][j-1] + d(j,i) < D[k][i] OR d(1,i) < D[k][i]
					else:
						if (j == 1):
							if d <= D[k][i]:
								D[k][i] = d
								B[k][i] = float(j)

						else:
							if (d+D[k-1][j-1] < D[k][i]):
								D[k][i] = d + D[k-1][j-1]
								B[k][i] = float(j)

	#Backtrack to find the clusters of the data points
	cluster_right = N
	cluster_left  = None
	Cluster 			= np.array([None]*N+1)   // record which cluster each point belongs to
	Centers 			= np.array([None]*K+1)	 // record the center of each cluster
	WithinSS 		 = np.array([None]*K+1)   // within sum of distance square of each cluster
	Size 				 = np.array([None]*K+1)	 // size of each cluster

	for k in reversed(range(1,K+1)):

		cluster_left = B[k][cluster_right]

		for i in range(cluster_left,cluster_right+1):
			Cluster[i] = k

		sum = 0.0
		for a in range(cluster_left,cluster_right+1):
			sum += x[a]

		Centers[k] = sum / (float(cluster_right) - cluster_left + 1)

		for a in range(cluster_left,cluster_right+1):
			WithinSS[k] += (x[a] - Centers[k]) **2

		Size[k] = cluster_right - int(cluster_left) + 1

		if (k > 1):
			cluster_right = int(cluster_left) - 1

	return {'Cluster':Cluster,'Centers':Centers,'Size':Size}
