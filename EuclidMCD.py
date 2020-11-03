import math
import argparse

DEBUG = None
ZERO_POLY = (0, )

# GCD(a, b) = d = x * a + y * b. Returns a tuple (d, x, y)
def euclid(a, b):
	u = (a, 1, 0)
	v = (b, 0, 1)

	q = u[0] // v[0]
	w = tuple(map(lambda x, y: x - y, u, map(lambda x: x * q, v)))

	while(w[0] != 0):
		if DEBUG:
			print('u = ' + str(u))
			print('v = ' + str(v))
			print('w = ' + str(w))
			print('')
		u = v
		v = w
		q = u[0] // v[0]
		w = tuple(map(lambda x, y: x - y, u, map(lambda x: x * q, v)))

	if DEBUG:
		print('u = ' + str(u))
		print('v = ' + str(v))
		print('w = ' + str(w))
		print('')
	return v

# Computes the inverse of x modulus modulus, if exists (i.e. GCD(x, modulus) = 1)
def mod_inverse(x, modulus):
	euclid_res = euclid_core(x, modulus)
	if euclid_res[0] != 1:
		raise ValueError(str(x) + " is not invertible modulus " + str(modulus))
	return mod_normalize(euclid_res[2], modulus)

# Normalize an integer number between 0 (included) and modulus (excluded)
def mod_normalize(x, modulus):
	orig = x

	if(x >= 0):
		x = x % modulus
	else:
		tmp = abs(x) // modulus
		if(abs(x) % modulus != 0):
			tmp = tmp + 1
		x = x + modulus * tmp

	#if DEBUG:
	#	print(str(orig) + " normalized to " + str(x) + " modulus " + str(modulus) + "\n")

	return x

def poly_string(f):
	ret = ''
	curr_deg = len(f) - 1
	for c in f:
		ret += str(c)
		if curr_deg > 1:
			ret += ' x^' + str(curr_deg)
		elif curr_deg == 1:
			ret += ' x'
		if curr_deg != 0:
			ret += ' + '
		curr_deg -= 1

	return ret

def remove_leading_zeros(f):
	i = 0
	while(f[i] == 0):
		i += 1
		if i >= len(f):
			return [0]

	return f[i:]


# Compute the polynomial subtraction in the ring (R[x], +, *)
def poly_sub(f, g):
	#if DEBUG:
	#	print("Computing the subtraction between " + poly_string(f) + " and " + poly_string(g) + "\n")
	lf = len(f)
	lg = len(g)
	if lf < lg:
		f = [0 for i in range(lg-lf)] + list(f)
	elif lg < lf:
		g = [0 for i in range(lf-lg)] + list(g)

	l = list(map(lambda x, y: 0 if math.isclose(x, y, abs_tol = 0.0001) else x - y, f, g))
	#if DEBUG:
	#	print('f = ' + poly_string(f))
	#	print('g = ' + poly_string(g))
	#	print('l = ' + poly_string(l))
	#	print('')

	return remove_leading_zeros(l)

# Compute polynomial division in the ring (R[x], +, *)
def poly_division(f, g):
	#if DEBUG:
	#	print("Computing division between " + poly_string(f) + " and " + poly_string(g) + "\n")
	deg_f = len(f) - 1
	deg_g = len(g) - 1
	res_deg = deg_f - deg_g
	q = [0 for x in range(res_deg + 1)]
	while(deg_f >= deg_g and f != [0]):
		sub = []
		curr = f[0] / g[0]
		curr_deg = deg_f - deg_g
		q[res_deg - curr_deg] = curr
		for i in range(deg_g + 1):
			sub.append(curr * g[i])
		for i in range(curr_deg):
			sub.append(0)
		
		#if DEBUG:
		#	print('curr = ' + str(curr))
		#	print('f = ' + poly_string(f))
		#	print('sub = ' + poly_string(sub))

		f = poly_sub(f, sub)
		deg_f = len(f) - 1

		#if DEBUG:
		#	print('new f = ' + poly_string(f))
		#	print('')
	if not q:
		q = [0]
	return (tuple(q), f)

def poly_mul(f, g):
	deg_f = len(f) - 1
	deg_g = len(g) - 1
	res_deg = deg_f + deg_g
	res = [0 for i in range(res_deg + 1)]

	iter_deg = deg_f if deg_f <= deg_g else deg_g
	iter_poly = f if iter_deg == deg_f else g
	other_deg = deg_g if iter_deg == deg_f else deg_f
	other_poly = g if iter_poly == f else f

	for i in range(iter_deg + 1):
		for j in range(other_deg + 1):
			local_mul = iter_poly[i] * other_poly[j]
			res[i + j] += local_mul

	return res

def poly_euclid(f, g):
	if DEBUG:
		print("Computing MCD(" + poly_string(f) + ", " + poly_string(g) + ")\n")
	u = (f, [1], [0])
	v = (g, [0], [1])

	q = poly_division(u[0], v[0])[0]

	w = tuple(map(lambda x, y: poly_sub(x, y), u, map(lambda x: poly_mul(x, q), v)))

	if DEBUG:
		print('u = ' + str(u))
		print('v = ' + str(v))
		print('w = ' + poly_string(w[0]))
		print('')

	while(w[0] != [0]):
		u = v
		v = w
		q = poly_division(u[0], v[0])[0]

		w = tuple(map(lambda x, y: poly_sub(x, y), u, map(lambda x: poly_mul(x, q), v)))

		if DEBUG:
			print('u = ' + str(u))
			print('v = ' + str(v))
			print('w = ' + poly_string(w[0]))
			print('')

	return v

def poly_mod_normalize(f, modulus):
	return tuple(map(lambda x: mod_normalize(x, modulus), f))

def poly_mod_division(f, g, modulus):
	#if DEBUG:
	#	print("Computing modular division between " + poly_string(f) + " and " + poly_string(g) + " with modulus " + str(modulus) + "\n")
	deg_f = len(f) - 1
	deg_g = len(g) - 1
	res_deg = deg_f - deg_g
	q = [0 for x in range(res_deg + 1)]
	while(deg_f >= deg_g and f != [0]):
		sub = []
		curr = f[0] * mod_inverse(g[0], modulus)
		curr_deg = deg_f - deg_g
		q[res_deg - curr_deg] = curr
		for i in range(deg_g + 1):
			sub.append(curr * g[i])
		for i in range(curr_deg):
			sub.append(0)
		
		#if DEBUG:
		#	print('curr = ' + str(curr))
		#	print('f = ' + poly_string(f))
		#	print('sub = ' + poly_string(sub))
		f = poly_mod_normalize(poly_sub(f, sub), modulus)
		f = remove_leading_zeros(f)
		#if DEBUG:
		#	print('new f = ' + poly_string(f))
		#	print('')
		deg_f = len(f) - 1

	if not q:
		q = [0]
	return (poly_mod_normalize(q, modulus), f)

def poly_euclid_core(f, g, modulus = None):
	if not modulus:
		return poly_euclid(f, g)

	f = poly_mod_normalize(f, modulus)
	g = poly_mod_normalize(g, modulus)

	u = (f, [1], [0])
	v = (g, [0], [1])

	q = poly_mod_division(u[0], v[0], modulus)[0]

	w = tuple(map(lambda x, y: poly_mod_normalize(poly_sub(x, y), modulus), u, map(lambda x: poly_mod_normalize(poly_mul(x, q), modulus), v)))

	if DEBUG:
		print('u = ' + str(u))
		print('v = ' + str(v))
		print('w = ' + poly_string(w[0]) + " " + str(w[1:]))
		print('')
	
	while(w[0] != ZERO_POLY):
		u = v
		v = w

		q = poly_mod_division(u[0], v[0], modulus)[0]
		w = tuple(map(lambda x, y: poly_mod_normalize(poly_sub(x, y), modulus), u, map(lambda x: poly_mod_normalize(poly_mul(x, q), modulus), v)))

		if DEBUG:
			print('u = ' + str(u))
			print('v = ' + str(v))
			print('w = ' + poly_string(w[0]) + " " + str(w[1:]))
			print('')

	return v


# If a modulus is passed, computations are modular, otherwise they are in R.
# Again, returns a tuple (d, x, y), where d = MCD(a, b) = x * a + y * b
def euclid_core(a, b, modulus = None):
	if not modulus:
		a = abs(a)
		b = abs(b)
		return euclid(a, b) if a > b else euclid(b, a)

	a = mod_normalize(a, modulus)
	b = mod_normalize(b, modulus)

	if b > a:
		tmp = b
		b = a
		a = tmp
		tmp = None

	u = (a, 1, 0)
	v = (b, 0, 1)

	q = u[0] * mod_inverse(v[0], modulus)
	w = tuple(map(lambda x, y: mod_normalize(x - y, modulus), u, map(lambda x: x * q, v)))

	while(w[0] != 0):
		if DEBUG:
			print('u = ' + str(u))
			print('v = ' + str(v))
			print('w = ' + str(w))
			print('')
		u = v
		v = w
		q = u[0] * mod_inverse(v[0], modulus)
		w = tuple(map(lambda x, y: mod_normalize(x - y, modulus), u, map(lambda x: x * q, v)))

	if DEBUG:
		print('u = ' + str(u))
		print('v = ' + str(v))
		print('w = ' + str(w))
		print('')
	return v

# Converts a string of type 'num / den' into the corresponding float num/den
def convert_to_float(s):
	arr = s.split('/')
	num = float(arr[0].strip())
	den = float(arr[1].strip())

	return num / den

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Euclid's algorithm to compute MCD")
	parser.add_argument('vals', nargs = 2)
	parser.add_argument('modulus', nargs = '?', type = int)
	parser.add_argument('-d', '--debug', action = 'store_true', default = False)
	parser.add_argument('--mod-inverse', action = 'store_true', default = False)
	parser.add_argument('--poly', action = 'store_true', default = False)
	parser.add_argument('--poly-div', action = 'store_true', default = False)
	parser.add_argument('--poly-mul', action = 'store_true', default = False)
	namespace = parser.parse_args()

	a = namespace.vals[0]
	b = namespace.vals[1]
	modulus = namespace.modulus
	DEBUG = namespace.debug

	if namespace.poly or namespace.poly_div or namespace.poly_mul:
		a = tuple(map(lambda s: convert_to_float(s.strip()) if '/' in s else float(s.strip()), a.split(',')))
		b = tuple(map(lambda s: convert_to_float(s.strip()) if '/' in s else float(s.strip()), b.split(',')))

		if namespace.poly:
			print(poly_string(poly_euclid_core(a, b, modulus)[0]))
		elif namespace.poly_div:
			div = poly_division(a, b) if not modulus else poly_mod_division(a, b, modulus)
			print("Quotient: " + poly_string(div[0]))
			print("Remainder: " + poly_string(div[1]))
		else:
			mul = poly_mul(a, b) if not modulus else poly_mod_normalize(poly_mul(a, b), modulus)
			print(poly_string(mul))
		exit(0)

	a = int(a)
	b = int(b)

	if namespace.mod_inverse:
		if(b < 0):
			print("Negative modulus, considering the absolute value")
			b = abs(b)
		print(mod_inverse(mod_normalize(a, b), b))
		exit(0)

	if(modulus and modulus < 0):
		print("Negative modulus, considering the absolute value")
		modulus = abs(int(modulus))
	
	x = euclid_core(a, b, modulus)
	print(x)