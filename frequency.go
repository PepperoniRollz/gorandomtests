package gorandomtests

import (
	"math"
)

func Frequency(n int) {
	var f, s_obs, p_value, sum, sqrt2 float64 = 0, 0, 0, 0, 1.41421356237309504880

	for i := 0; i < n; i++ {
		sum += 2*epsilon[i] - 1
	}
	s_obs = math.Abs(sum) / math.Sqrt(float64(n))
	f = s_obs / sqrt2
	p_value = math.Erfc(f)

}
