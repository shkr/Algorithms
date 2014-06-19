#in KCenter1D.jl

function KCenter1D(y::Array{Float64,1},K::Int64)

	# Description
	# -----------
	# one-dimensional cluster algorithm implemented in julia
	# y is input one-dimensional vector and K stands for the cluster level

	# Input
	# -----------
	# y -- a vector of sorted numbers
	# K -- the number of clusters expected

	#Trivial case
	# -----------
	# K > |set(y)|
	# |set(y)| == 1
	if (K>length(unique(y)) || length(unique(y))==1)
		return (ASCIIString=>Array{Float64})["Cluster"=>ones(Int32,length(y)),"Centers"=>y,"Size"=>ones(Int32,K),"WithinSS"=>zeros(Float64,length(y))]
	end

	# Assumptions
	# -----------
	# K < |set(y)|
	# |set(y)| > 1


	N = length(y) # N: is the size of input vector
	

	#D - Matrix rows 1 to K, columns 1 to N, OPTIMAL Solution WithinSS for KCenter(y,K)
	D = Array(Float64,(K,N))
	#B - Matrix rows 1 to K, columns 1 to N, Position of leftmost element of Kth cluster of OPTIMAL Solution
	B = Array(Float64,(K,N))

	#Base case : Set for all K rows, Column 1 to 0
	B[1:K,1] = 1.0
	D[1:K,1] = 0.0
	
	# running values
	# μ_y1     :     mean of when k=1
	# μ_yj	   :     mean of rightmost-cluster when k>1
	# d        :     SS (Sum of Squares distance) from mean for the rightmost-cluster

	for k=1:K

		μ_y1 = y[1]

		for i=2:N

			# when k = 1, D(1,i) = d(1 to i)
			if (k == 1)

				D[1,i] = D[1,i-1] + ((float64(i-1)/float64(i))*(y[i]-μ_y1)^2)
				μ_y1   = (float64(i-1)*μ_y1 + y[i])/float64(i)
				B[1,i] = 1

			#when k > 1, D(k,i)  = min ( D(k-1,j-1) + d(j,i) ) for j>=1
			else

				d      = 0.0   			#Initiallly no item in rightmost-cluster
				μ_yj   = 0.0	 		#Initiallly no item in rightmost-cluster
				D[k,i] = -1  			#value of D(k,i) is not set

				for j in i:-1:1

					#One-by-one add item (j,i) to rightmost-cluster
					d += (float64(i - j) / float64(i-j+1)) * (y[j] - μ_yj)^2   #update d(j,i)
					μ_yj = (y[j] + float64(i-j)*μ_yj) / float64(i-j+1)	   	   #update mean(j,i)

					#initialize D(k,i) with D(k,i) = D(k-1,i-1) + d
					if (D[k,i] == -1)
						if (j == 1)
							D[k,i] = d
							B[k,i] = float64(j)
						else
							D[k,i] = d + D[k-1,j-1]
							B[k,i] = float64(j)
						end
					#If already initialized find if D[k-1][j-1] + d(j,i) < D[k][i] OR d(1,i) < D[k][i]
					else
						if (j == 1)
							if d <= D[k,i]
								D[k,i] = d
								B[k,i] = float64(j)
							end
						else
							if (d+D[k-1,j-1] < D[k,i])
								D[k,i] = d + D[k-1,j-1]
								B[k,i] = float64(j)
							end
						end

					end

				end

			end

		end

	end

	#Backtrack to find the clusters of the data points
	cluster_right = N
	cluster_left  = -1
	Cluster 			= Array(Int32,N)     # record which cluster each point belongs to
	Center   			= Array(Float64,K)	 # record the center of each cluster
	WithinSS 		 	= Array(Float64,K)   # within sum of distance square of each cluster
	Size 				= Array(Float64,K)	 # size of each cluster

	for k=K:-1:1

		cluster_left = B[k,cluster_right]

		for i=cluster_left:cluster_right
			Cluster[i] = k
		end

		ss = 0.0
		
		for a=cluster_left:cluster_right
			ss += y[a]
		end

		Center[k] = ss / (float64(cluster_right) - cluster_left + 1)

		for a=cluster_left:cluster_right
			WithinSS[k] += (y[a] - Center[k])^2
		end

		Size[k] = cluster_right - int(cluster_left) + 1

		if (k > 1)
			cluster_right = int(cluster_left) - 1
		end
	
	end

	return (ASCIIString=>Array{Float64})["Cluster"=>Cluster,"Centers"=>Center,"Size"=>Size,"WithinSS"=>WithinSS]
end