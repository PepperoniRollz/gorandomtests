package gorandomtests

func BlockFrequency(M, n int) {

	var i, j, N, blockSum int
	var p_value, sum, pi, v, chi_squared float64

	n = n / M
	sum = 0.0

	for i := 0; i < N; i++ {
		blockSum = 0
		for j := 0; j < M; j++ {
			blockSum += epsilon[j+i*M]
		}
		pi = float64(blockSum) / float64(M)
		v = pi - 0.5
		sum += v * v
	}

	chi_squared = 4.0 * M * sum
	p_value = cephes_igamc(N/2.0, chi_squared/2.0)

}
