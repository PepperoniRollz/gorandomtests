package gorandomtests

import (
	"fmt"
	"math"
)

// Global variables
var (
	relError = 1e-12
	MACHEP   = 1.11022302462515654042e-16  // 2**-53
	MAXLOG   = 7.09782712893383996732224e2 // log(MAXNUM)
	MAXNUM   = 1.7976931348623158e308      // 2**1024*(1-MACHEP)
	PI       = 3.14159265358979323846      // pi, duh!
	big      = 4.503599627370496e15
	biginv   = 2.22044604925031308085e-16
	sgngam   = 0
)

// cephes_igamc computes the complementary incomplete gamma function.
func cephesIgamc(a, x float64) float64 {

	var ans, ax, c, yc, r, t, y, z float64
	var pk, pkm1, pkm2, qk, qkm1, qkm2 float64

	if (x <= 0) || (a <= 0) {
		return (1.0)
	}

	if (x < 1.0) || (x < a) {
		return (1.e0 - cephesIgam(a, x))
	}

	ax = a*math.Log(x) - x - cephesLgam(a)

	if ax < -MAXLOG {
		fmt.Println("igamc: UNDERFLOW")
		return 0.0
	}
	ax = math.Exp(ax)

	/* continued fraction */
	y = 1.0 - a
	z = x + y + 1.0
	c = 0.0
	pkm2 = 1.0
	qkm2 = x
	pkm1 = x + 1.0
	qkm1 = z * x
	ans = pkm1 / qkm1

	for {
		c += 1.0
		y += 1.0
		z += 2.0
		yc = y * c
		pk = pkm1*z - pkm2*yc
		qk = qkm1*z - qkm2*yc
		if qk != 0 {
			r = pk / qk
			t = math.Abs((ans - r) / r)
			ans = r
		} else {
			t = 1.0
		}
		pkm2, pkm1 = pkm1, pk
		qkm2, qkm1 = qkm1, qk
		if math.Abs(pk) > big {
			pkm2 *= biginv
			pkm1 *= biginv
			qkm2 *= biginv
			qkm1 *= biginv
		}
		if t <= MACHEP {
			break // Exit the loop when the condition is false
		}
	}

	return ans * ax

}

func cephesIgam(a, x float64) float64 {
	if x <= 0 || a <= 0 {
		return 0.0
	}

	if x > 1.0 && x > a {
		return 1.0 - cephesIgamc(a, x)
	}

	ax := a*math.Log(x) - x - cephesLgam(a)
	if ax < -MAXLOG {
		fmt.Println("igam: UNDERFLOW")
		return 0.0
	}
	ax = math.Exp(ax)

	// Power series
	r := a
	c := 1.0
	ans := 1.0

	for c/ans > MACHEP {
		r += 1.0
		c *= x / r
		ans += c
	}

	return ans * ax / a
}

var A = []uint16{
	0x6661, 0x2733, 0x9850, 0x3f4a,
	0xe943, 0xb580, 0x7fbd, 0xbf43,
	0x5ebb, 0x20dc, 0x019f, 0x3f4a,
	0xa5a1, 0x16b0, 0xc16c, 0xbf66,
	0x554b, 0x5555, 0x5555, 0x3fb5,
}

var B = []uint16{
	0x6761, 0x8ff3, 0x8901, 0xc095,
	0xb93e, 0x355b, 0xf234, 0xc0e2,
	0x89e5, 0xf890, 0x3d73, 0xc114,
	0xdb51, 0xf994, 0xbc82, 0xc131,
	0xf20b, 0x0219, 0x4589, 0xc13a,
	0x055e, 0x5418, 0x0c67, 0xc12a,
}

var C = []uint16{
	// 0x0000, 0x0000, 0x0000, 0x3ff0,
	0x12b2, 0x1cf3, 0xfd0d, 0xc075,
	0xd757, 0x7b89, 0xaa0d, 0xc0d0,
	0x4c9b, 0xb974, 0xeb84, 0xc10a,
	0x0043, 0x7195, 0x6286, 0xc131,
	0xf34c, 0x892f, 0x5255, 0xc143,
	0xe14a, 0x6a11, 0xce4b, 0xc13e,
}

const MAXLGM = 2.556348e305

func cephesLgam(x float64) (result float64) {
	var p, q, u, w, z float64

	sgngam = 1

	if x < -34.0 {
		q = -x
		w = cephesLgam(q) // Recursive call modifies sgngam!
		p = math.Floor(q)
		if p == q {
			fmt.Println("lgam: Singularity")
			return float64(sgngam) * MAXNUM
		}
		i := int(p)
		sgngam = 1 - 2*(i&1) // If i is even, sgngam = -1; else sgngam = 1

		z = q - p
		if z > 0.5 {
			p += 1.0
			z = p - q
		}
		z = q * math.Sin(PI*z)
		if z == 0.0 {
			fmt.Println("lgam: Singularity")
			return float64(sgngam) * MAXNUM
		}
		return math.Log(PI) - math.Log(z) - w
	}

	if x < 13.0 {
		z = 1.0
		p = 0.0
		u = x
		for u >= 3.0 {
			p -= 1.0
			u = x + p
			z *= u
		}
		for u < 2.0 {
			if u == 0.0 {
				fmt.Println("lgam: Singularity")
				return float64(sgngam) * MAXNUM
			}
			z /= u
			p += 1.0
			u = x + p
		}
		if z < 0.0 {
			sgngam = -1
			z = -z
		} else {
			sgngam = 1
		}
		if u == 2.0 {
			return math.Log(z)
		}
		p -= 2.0
		x += p
		p = x * cephesPolevl(x, B, 5) / cephesP1evl(x, C, 6)
		result = math.Log(z) + p
		// Here you would include the actual polynomial evaluation logic
		return result // Update this return based on your polynomial evaluation
	}

	if x > MAXLGM {
		fmt.Println("lgam: OVERFLOW")
		return float64(sgngam) * MAXNUM
	}

	q = (x-0.5)*math.Log(x) - x + math.Log(math.Sqrt(2*PI))
	if x <= 1.0e8 {
		p = 1.0 / (x * x)
		if x >= 1000.0 {
			q += ((7.9365079365079365079365e-4*p-2.7777777777777777777778e-3)*p + 0.0833333333333333333333) / x
		} else {
			q += cephesPolevl(p, A, 4) / x
		}
	}

	return q
}

// cephesPolevl evaluates a polynomial given coefficients and a value for x.
func cephesPolevl(x float64, coef []uint16, N int) float64 {
	ans := float64(coef[0])
	for i := 1; i <= N; i++ {
		ans = ans*x + float64(coef[i])
	}
	return ans
}

// cephesP1evl evaluates a polynomial given coefficients and a value for x, skipping the first term.
func cephesP1evl(x float64, coef []uint16, N int) float64 {
	ans := x + float64(coef[0])
	for i := 1; i < N; i++ {
		ans = ans*x + float64(coef[i])
	}
	return ans
}

// cephesErf computes the error function of x.
func cephesErf(x float64) float64 {
	const twoSqrtPi = 1.128379167095512574
	if math.Abs(x) > 2.2 {
		return 1.0 - cephesErfc(x)
	}
	sum, term, xsqr := x, x, x*x
	for j := 1; math.Abs(term)/sum > relError; j += 2 {
		term *= xsqr / float64(j)
		sum -= term / float64(2*j+1)
		term *= xsqr / float64(j+1)
		sum += term / float64(2*j+3)
	}
	return twoSqrtPi * sum
}

// cephesErfc computes the complementary error function of x.
func cephesErfc(x float64) float64 {
	const oneSqrtPi = 0.564189583547756287
	if math.Abs(x) < 2.2 {
		return 1.0 - cephesErf(x)
	}
	if x < 0 {
		return 2.0 - cephesErfc(-x)
	}
	a, b, c, d := 1.0, x, x, x*x+0.5
	q2 := b / d
	for n, q1 := 1.0, 0.0; math.Abs(q1-q2)/q2 > relError; n += 0.5 {
		t := a*n + b*x
		a, b = b, t
		t = c*n + d*x
		c, d = d, t
		q1 = q2
		q2 = b / d
	}
	return oneSqrtPi * math.Exp(-x*x) * q2
}

// cephesNormal computes the normal distribution function for x.
func cephesNormal(x float64) float64 {
	const sqrt2 = 1.414213562373095048801688724209698078569672
	var arg, result float64
	if x > 0 {
		arg = x / sqrt2
		result = 0.5 * (1 + cephesErf(arg))
	} else {
		arg = -x / sqrt2
		result = 0.5 * (1 - cephesErf(arg))
	}
	return result
}
