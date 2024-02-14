package gorandomtests

import (
	"math"
)

func runs(n int) {
	var S int
	var pi, V, erfc_arg, p_value float64

	S = 0

	for k := 0; k < n; k++ {
		if epsilon[k] != 0 {
			S++
		}
	}
	pi = float64(S) / float64(n)

	if math.Abs(pi-0.5) > (2.0 / math.Sqrt(n)) {
		p_value = 0.0
	} else {
		V = 1
		for k := 1; k < n; k++ {
			if epsilon[k] != epsilon[k-1] {
				V++
			}
		}

		erfc_arg = math.Abs(V-2.0*float64(n)*pi*(1-pi)) / (2.0 * pi * (1 - pi) * math.Sqrt(2*float64(n)))
		p_value = math.Erfc(erfc_arg)
	}

}

// int		S, k;
// double	pi, V, erfc_arg, p_value;

// S = 0;
// for ( k=0; k<n; k++ )
// 	if ( epsilon[k] )
// 		S++;
// pi = (double)S / (double)n;

// if ( fabs(pi - 0.5) > (2.0 / sqrt(n)) ) {
// 	fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
// 	fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
// 	fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
// 	p_value = 0.0;
// }
// else {

// 	V = 1;
// 	for ( k=1; k<n; k++ )
// 		if ( epsilon[k] != epsilon[k-1] )
// 			V++;

// 	erfc_arg = fabs(V - 2.0 * n * pi * (1-pi)) / (2.0 * pi * (1-pi) * sqrt(2*n));
// 	p_value = erfc(erfc_arg);
